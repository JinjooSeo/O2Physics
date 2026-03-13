// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "PWGDQ/Core/MixingHandler.h"
#include "PWGDQ/Core/MixingLibrary.h"
#include "PWGDQ/Core/VarManager.h"

#include <TObjArray.h>

#include "rapidjson/document.h"
#include "rapidjson/error/en.h"

#include <cstring>
#include <iostream>
#include <fstream>
#include <memory>
using namespace std;

#include <TMath.h>
#include <TTimeStamp.h>
#include <TRandom.h>

ClassImp(MixingHandler);

namespace
{
int getVariableFromName(TString variableName)
{
  variableName = variableName.Strip(TString::kBoth, ' ');
  if (variableName.IsNull()) {
    return VarManager::kNothing;
  }

  if (VarManager::fgVarNamesMap.empty()) {
    VarManager::SetDefaultVarNames();
  }

  auto variable = VarManager::fgVarNamesMap.find(variableName);
  if (variable != VarManager::fgVarNamesMap.end()) {
    return variable->second;
  }

  if (!variableName.BeginsWith("k")) {
    TString prefixedName = Form("k%s", variableName.Data());
    variable = VarManager::fgVarNamesMap.find(prefixedName);
    if (variable != VarManager::fgVarNamesMap.end()) {
      return variable->second;
    }
  }

  return VarManager::kNothing;
}

int inferVariableFromMixingName(TString mixingName)
{
  mixingName = mixingName.Strip(TString::kBoth, ' ');
  if (mixingName.BeginsWith("CentralityFT0C")) {
    return VarManager::kCentFT0C;
  }
  if (mixingName.BeginsWith("Centrality")) {
    return VarManager::kCentVZERO;
  }
  if (mixingName.BeginsWith("Mult")) {
    return VarManager::kVtxNcontrib;
  }
  if (mixingName.BeginsWith("Vtx")) {
    return VarManager::kVtxZ;
  }
  if (mixingName.BeginsWith("Occupancy")) {
    return VarManager::kTrackOccupancyInTimeRange;
  }
  if (mixingName.BeginsWith("Psi2A")) {
    return VarManager::kPsi2A;
  }
  if (mixingName.BeginsWith("Psi2B")) {
    return VarManager::kPsi2B;
  }
  if (mixingName.BeginsWith("Psi2C")) {
    return VarManager::kPsi2C;
  }
  if (mixingName.BeginsWith("MedianTimeA")) {
    return VarManager::kNTPCmedianTimeLongA;
  }
  if (mixingName.BeginsWith("PileUpA")) {
    return VarManager::kNTPCcontribLongA;
  }
  return VarManager::kNothing;
}

bool validateBinLimits(const rapidjson::Value& binLimits, const char* mixingName)
{
  if (!binLimits.IsArray()) {
    LOG(fatal) << "Mixing definition " << mixingName << " must provide bin limits as an array";
    return false;
  }
  if (binLimits.GetArray().Size() < 2) {
    LOG(fatal) << "Mixing definition " << mixingName << " must provide at least two bin edges";
    return false;
  }

  double previousEdge = 0.0;
  bool firstEdge = true;
  for (const auto& edge : binLimits.GetArray()) {
    if (!edge.IsNumber()) {
      LOG(fatal) << "Mixing definition " << mixingName << " contains a non-numeric bin edge";
      return false;
    }
    double currentEdge = edge.GetDouble();
    if (!firstEdge && currentEdge <= previousEdge) {
      LOG(fatal) << "Mixing definition " << mixingName << " must provide strictly increasing bin edges";
      return false;
    }
    previousEdge = currentEdge;
    firstEdge = false;
  }

  return true;
}

int getJSONMixingVariable(const rapidjson::Value& mixing, const char* mixingName)
{
  if (mixing.IsObject() && mixing.HasMember("var")) {
    const auto& varField = mixing.FindMember("var")->value;
    if (!varField.IsString()) {
      LOG(fatal) << "Mixing definition " << mixingName << " has a non-string var field";
      return VarManager::kNothing;
    }
    int variable = getVariableFromName(varField.GetString());
    if (variable == VarManager::kNothing) {
      LOG(fatal) << "Mixing definition " << mixingName << " uses an unknown variable " << varField.GetString();
      return VarManager::kNothing;
    }
    return variable;
  }

  int inferredVariable = inferVariableFromMixingName(mixingName);
  if (inferredVariable == VarManager::kNothing) {
    LOG(fatal) << "Mixing definition " << mixingName << " must specify a valid var field";
  }
  return inferredVariable;
}

bool validateJSONMixingDefinition(const rapidjson::Value& mixing, const char* mixingName)
{
  if (mixing.IsArray()) {
    if (inferVariableFromMixingName(mixingName) == VarManager::kNothing) {
      LOG(fatal) << "Cannot infer the variable for mixing definition " << mixingName << "; please specify var explicitly";
      return false;
    }
    return validateBinLimits(mixing, mixingName);
  }

  if (!mixing.IsObject()) {
    LOG(fatal) << "Mixing definition " << mixingName << " must be an object or an array of bin edges";
    return false;
  }

  if (!mixing.HasMember("binLimits")) {
    LOG(fatal) << "Mixing definition " << mixingName << " is missing the binLimits field";
    return false;
  }

  (void)getJSONMixingVariable(mixing, mixingName);
  return validateBinLimits(mixing.FindMember("binLimits")->value, mixingName);
}

bool addJSONMixingVariable(MixingHandler* mh, const rapidjson::Value& mixing, const char* mixingName)
{
  if (!validateJSONMixingDefinition(mixing, mixingName)) {
    return false;
  }

  int variable = getJSONMixingVariable(mixing, mixingName);
  const auto& binLimitsJSON = mixing.IsArray() ? mixing : mixing.FindMember("binLimits")->value;
  std::vector<float> binLimits;
  binLimits.reserve(binLimitsJSON.GetArray().Size());
  for (const auto& edge : binLimitsJSON.GetArray()) {
    binLimits.push_back(static_cast<float>(edge.GetDouble()));
  }

  mh->AddMixingVariable(variable, binLimits.size(), binLimits);
  return true;
}
} // namespace

//_________________________________________________________________________
MixingHandler::MixingHandler() : TNamed(),
                                 fIsInitialized(kFALSE),
                                 fVariableLimits(),
                                 fVariables()
{
  //
  // default constructor
  //
}

//_________________________________________________________________________
MixingHandler::MixingHandler(const char* name, const char* title) : TNamed(name, title),
                                                                    fIsInitialized(kFALSE),
                                                                    fVariableLimits(),
                                                                    fVariables()
{
  //
  // Named constructor
  //
}

//_________________________________________________________________________
MixingHandler::~MixingHandler()
{
  //
  // destructor
  //
}

//_________________________________________________________________________
void MixingHandler::AddMixingVariable(int var, int nBins, float* binLims)
{
  //
  // add a mixing variable
  //
  fVariables.push_back(var);
  TArrayF varBins;
  varBins.Set(nBins, binLims);
  fVariableLimits.push_back(varBins);
  VarManager::SetUseVariable(var);
}

//_________________________________________________________________________
void MixingHandler::AddMixingVariable(int var, int nBins, std::vector<float> binLims)
{

  float* bins = new float[nBins];
  for (int i = 0; i < nBins; ++i) {
    bins[i] = binLims[i];
  }
  AddMixingVariable(var, nBins, bins);
}

//_________________________________________________________________________
int MixingHandler::GetMixingVariable(VarManager::Variables var)
{
  int i = 0;
  for (auto v = fVariables.begin(); v != fVariables.end(); v++, i++) {
    if (*v == var) {
      return i;
    }
  }
  return -1;
}

//_________________________________________________________________________
std::vector<float> MixingHandler::GetMixingVariableLimits(VarManager::Variables var)
{
  std::vector<float> binLimits;
  int i = 0;
  for (auto v = fVariables.begin(); v != fVariables.end(); v++, i++) {
    if (*v == var) {
      for (int iBin = 0; iBin < fVariableLimits[i].GetSize(); ++iBin) {
        binLimits.push_back(fVariableLimits[i].At(iBin));
      }
      break;
    }
  }
  return binLimits;
}

//_________________________________________________________________________
void MixingHandler::Init()
{
  //
  // Initialization of pools
  //       The correct event category will be retrieved using the function FindEventCategory()
  //
  int size = 1;
  for (auto v : fVariableLimits) {
    size *= (v.GetSize() - 1);
  }
  (void)size;
  fIsInitialized = kTRUE;
}

//_________________________________________________________________________
int MixingHandler::FindEventCategory(float* values)
{
  //
  // Find the event category corresponding to the added mixing variables
  //
  if (fVariables.size() == 0) {
    return -1;
  }
  if (!fIsInitialized) {
    Init();
  }

  std::vector<int> bin;
  int iVar = 0;
  for (auto v = fVariableLimits.begin(); v != fVariableLimits.end(); v++, iVar++) {
    int binValue = TMath::BinarySearch((*v).GetSize(), (*v).GetArray(), values[fVariables[iVar]]);
    bin.push_back(binValue);
    if (bin[iVar] == -1 || bin[iVar] == (*v).GetSize() - 1) {
      return -1; // all variables must be inside limits
    }
  }

  int category = 0;
  int tempCategory = 1;
  int iv1 = 0;
  int iv2 = 0;
  for (auto v1 = fVariables.begin(); v1 != fVariables.end(); v1++, iv1++) {
    tempCategory = 1;
    iv2 = iv1;
    for (auto v2 = v1; v2 != fVariables.end(); v2++, iv2++) {
      if (iv2 == iv1) {
        tempCategory *= bin[iv2];
      } else {
        tempCategory *= (fVariableLimits[iv2].GetSize() - 1);
      }
    }
    category += tempCategory;
  }
  return category;
}

//_________________________________________________________________________
int MixingHandler::GetBinFromCategory(VarManager::Variables var, int category) const
{
  //
  // find the bin in variable var for the n-dimensional "category"
  //
  if (fVariables.size() == 0) {
    return -1;
  }

  // Search for the position of the variable "var" in the internal variable list of the handler
  int tempVar = 0;
  for (auto v = fVariables.begin(); v != fVariables.end(); v++, tempVar++) {
    if (*v == var) {
      break;
    }
  }

  // extract the bin position in variable "var" from the category
  int norm = 1;
  for (int i = fVariables.size() - 1; i > tempVar; --i) {
    norm *= (fVariableLimits[i].GetSize() - 1);
  }
  int truncatedCategory = category - (category % norm);
  truncatedCategory /= norm;
  return truncatedCategory % (fVariableLimits[tempVar].GetSize() - 1);
}

void o2::aod::dqmixing::AddMixingVariables(MixingHandler* mh, const char* mixingVariables, const char* json)
{
  if (!mh) {
    LOG(fatal) << "MixingHandler pointer is null";
    return;
  }

  TString mixVarsString = mixingVariables ? mixingVariables : "";
  mixVarsString = mixVarsString.Strip(TString::kBoth, ' ');
  if (mixVarsString.Length() == 0) {
    return;
  }

  rapidjson::Document document;
  bool hasJSONMixingDefinitions = false;
  TString jsonString = json ? json : "";
  jsonString = jsonString.Strip(TString::kBoth, ' ');
  if (jsonString.Length() > 0) {
    LOG(info) << "========================================== interpreting JSON for mixing variables";
    LOG(info) << "      json string is: " << json;

    rapidjson::ParseResult ok = document.Parse(json);
    if (!ok) {
      LOG(fatal) << "JSON parse error: " << rapidjson::GetParseErrorFunc(ok.Code()) << " (" << ok.Offset() << ")";
      TString str = "";
      for (int i = ok.Offset() - 30; i < static_cast<int>(ok.Offset()) + 50; i++) {
        if ((i >= 0) && (i < static_cast<int>(strlen(json)))) {
          str += json[i];
        }
      }
      LOG(fatal) << "**** Parsing error is somewhere here: " << str.Data();
      return;
    }
    if (!document.IsObject()) {
      LOG(fatal) << "Mixing JSON must be a top-level object keyed by names used in cfgMixingVars";
      return;
    }
    hasJSONMixingDefinitions = true;
  }

  std::unique_ptr<TObjArray> objArray(mixVarsString.Tokenize(","));
  if (!objArray) {
    return;
  }

  for (int iVar = 0; iVar < objArray->GetEntries(); ++iVar) {
    TString mixingName = objArray->At(iVar)->GetName();
    mixingName = mixingName.Strip(TString::kBoth, ' ');
    if (mixingName.Length() == 0) {
      continue;
    }

    int nMixingVariablesBefore = mh->GetNMixingVariables();
    SetUpMixing(mh, mixingName.Data());
    if (mh->GetNMixingVariables() > nMixingVariablesBefore) {
      continue;
    }

    if (hasJSONMixingDefinitions) {
      auto jsonMixing = document.FindMember(mixingName.Data());
      if (jsonMixing != document.MemberEnd()) {
        if (addJSONMixingVariable(mh, jsonMixing->value, mixingName.Data())) {
          LOG(info) << "Configured mixing variable " << mixingName.Data() << " from JSON";
          continue;
        }
      }
    }

    LOG(fatal) << "Did not find mixing configuration " << mixingName.Data() << " in MixingLibrary or cfgMixingVarsJSON";
  }
}

#include <XGBoostPredictor.h>

#include <limits>
#include <sstream>


XGBoostPredictorException::XGBoostPredictorException(
    std::string const &functionName, int exitCode) {
  std::ostringstream msg;
  msg << "Call to XGBoost C API function \"" << functionName
      << "\" terminated with error code " << exitCode << ".";
  message_ = msg.str();
}


XGBoostPredictor::XGBoostPredictor(std::string const &path, int numFeatures)
    : booster_{nullptr}, numFeatures_{numFeatures} {
  CheckCall("XGBoosterCreate", XGBoosterCreate(nullptr, 0, &booster_));
  CheckCall("XGBoosterSetParam", XGBoosterSetParam(booster_, "nthread", "1"));
  CheckCall("XGBoosterLoadModel", XGBoosterLoadModel(booster_, path.c_str()));
}


XGBoostPredictor::~XGBoostPredictor() {
  XGBoosterFree(booster_);
}


float XGBoostPredictor::Predict(float const *x) const {
  DMatrixHandle dmat{nullptr};
  CheckCall(
      "XGDMatrixCreateFromMat",
      XGDMatrixCreateFromMat(
          x, 1, numFeatures_, std::numeric_limits<float>::quiet_NaN(), &dmat));

  bst_ulong numScores = 0;
  float const *scores{nullptr};
  CheckCall(
      "XGBoosterPredict",
      XGBoosterPredict(booster_, dmat, 0, 0, &numScores, &scores));
  float const prediction = scores[0];

  CheckCall("XGDMatrixFree", XGDMatrixFree(dmat));
  return prediction;
}


void XGBoostPredictor::CheckCall(
    std::string const &functionName, int exitCode) {
  if (exitCode != 0)
    throw XGBoostPredictorException(functionName, exitCode);
}


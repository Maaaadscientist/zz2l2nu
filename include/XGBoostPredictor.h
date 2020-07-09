#ifndef HZZ2L2NU_INCLUDE_XGBOOSTPREDICTOR_H_
#define HZZ2L2NU_INCLUDE_XGBOOSTPREDICTOR_H_

#include <exception>
#include <string>

#include <xgboost/c_api.h>


/// Exception class for calls to XGBoost C API
class XGBoostPredictorException : public std::exception {
 public:
  /**
   * \brief Constructor
   *
   * \param[in] functionName  Name of the C API function that was called.
   * \param[in] exitCode  Exit code returned by the function.
   */
  XGBoostPredictorException(std::string const &functionName, int exitCode);

  /// Returns the explanatory string
  char const *what() const noexcept override {
    return message_.c_str();
  }

 private:
  /// Explanatory string
  std::string message_;
};


/**
 * \brief Simple wrapper around XGBoost C API
 *
 * This class loads a saved XGBoost model from a file and computes its
 * predictions, for one example at a time. The predictions require memory
 * allocation and deallocation for each example, which is not optimal.
 * Unfortunately, in XGBoost 0.90 API this is the only way to compute
 * predictions per example.
 */
class XGBoostPredictor {
 public:
  /**
   * \brief Constructor
   *
   * \param[in] path  Path to a file with XGBoost model.
   * \param[in] numFeatures  Number of input features.
   */
  XGBoostPredictor(std::string const &path, int numFeatures);

  ~XGBoostPredictor();

  /**
   * \brief Computes prediction of the XGBoost model for a single example
   *
   * \param[in] x  Array of features describing the example. Starting from the
   *   given position, \c numFeatures_ elements will be read.
   * \return Prediction of the model.
   */
  float Predict(float const *x) const;

 private:
  /**
   * \brief Throws XGBoostPredictorException if the exit code is not 0
   *
   * Arguments are forwarded to XGBoostPredictorException.
   */
  static void CheckCall(std::string const &functionName, int exitCode);

  /// Pointer to XGBoost model
  BoosterHandle booster_;

  /// Number of input features specified in the constructor
  int numFeatures_;
};

#endif  // HZZ2L2NU_INCLUDE_XGBOOSTPREDICTOR_H_


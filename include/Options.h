#ifndef OPTIONS_H_
#define OPTIONS_H_

#include <initializer_list>
#include <sstream>
#include <stdexcept>
#include <string>

#include <boost/program_options.hpp>
#include <yaml-cpp/yaml.h>

#include <Logger.h>


/**
 * \brief Provides access to command line options and configuration file
 *
 * Implements the parsing with Boost.Program_options. Values of options can be
 * accessed by their labels. All possible options must be registered beforehand.
 *
 * The configuration file must be of YAML format. It is parsed with
 * <a href="https://github.com/jbeder/yaml-cpp">yaml-cpp</library>.
 */
class Options {
 public:
  /**
   * \brief Exception for class Options
   */
  class Error : public std::runtime_error {
   public:
    Error(std::string const &message) : std::runtime_error{message} {};
  };

  using Group = boost::program_options::options_description;
  
  /**
   * \brief Constructor from command line arguments
   *
   * \param[in] argc,argv  Number of command line arguments and an array with
   *   their values, as given to \c main.
   * \param[in] optionGroups  Zero or more instances of
   *   boost::program_options::options_description that describe possible
   *   command line options.
   *
   * If an unregistered option is encountered, terminates the program. Several
   * options are added automatically:
   *  - \c -h,--help  Prints usage information and exists the program.
   *  - \c --config   Sets location of the configuration file. The path is
   *    resolved using FileInPath.
   *  - \c --version  Prints version and exits the program.
   *  - \c -v,--verbosity  Sets the verbosity level for the log.
   */
  Options(int argc, char **argv,
          std::initializer_list<Group> const &optionGroups);

  /**
   * \brief Checks if an option with the given label exists
   */
  bool Exists(std::string const &label) const;

  /**
   * \brief Returns the value of the option with the given label
   *
   * The value is represented with the given type. If the requested option has
   * not been specified in the call to the program or if the option does not
   * have a value (i.e. this is a pure flag), throws an exception of type
   * Options::Error.
   */
  template<typename T>
  T GetAs(std::string const &label) const;

  /**
   * \brief Checks and returns the value of the option with the given label
   *
   * \param[in] label    Label that identifies the option.
   * \param[in] checker  Object of type Checker that defines
   *   <tt>bool operator()(T const &)</tt> and checks the validity of the
   *   returned value for the option.
   * 
   * Behaves in the same way as \ref GetAs but additionally checks the value of
   * the option using the provided functor. If the value fails the check, throws
   * an exception of type Options::Error. In a typical use case the checker will
   * be a lambda function.
   */
  template<typename T, typename Checker>
  T GetAsChecked(std::string const &label, Checker const &checker) const;

  /**
   * \brief Returns parsed YAML configuration
   *
   * A brief introduction to the parsing library is available
   * <a href="https://github.com/jbeder/yaml-cpp/wiki/Tutorial">here</a>.
   * If no configuration file has been given, this method throws an exception of
   * type Options::Error.
   */
  YAML::Node const &GetConfig() const;

 private:
  /**
   * \brief Prints usage instructions
   */
  void PrintUsage() const;

  /// Name of the program extracted from the first command line argument
  std::string programName_;

  /// All registered options
  boost::program_options::options_description allOptions_;

  /// Map with parsed options
  boost::program_options::variables_map optionMap_;

  /**
   * \brief Parsed YAML configuration
   *
   * If no configuration file has been given in command line arguments, the
   * type of this node is Null.
   */
  YAML::Node config_;
};


template<typename T>
T Options::GetAs(std::string const &label) const {
  if (not Exists(label)) {
    std::ostringstream message;
    message << "Unknown option \"" << label << "\"";
    LOG_ERROR << message.str();
    throw Error(message.str());
  }

  return optionMap_[label].as<T>();
}


template<typename T, typename Checker>
T Options::GetAsChecked(std::string const &label,
                        Checker const &checker) const {
  T const value = GetAs<T>(label);

  if (not checker(value)) {
    std::ostringstream message;
    message << "Invalid value read for option \"" << label << "\": " << value;
    LOG_ERROR << message.str();
    throw Error(message.str());
  }

  return value;
}

#endif  // OPTIONS_H_



#include "catch.hpp"

namespace Catch
{
  struct MPIReporter : ConsoleReporter {

    MPIReporter(ReporterConfig const &_config) : ConsoleReporter(_config) {}

    virtual ~MPIReporter();

    static std::string getDescription() {
      return "Reports test results for MPI unit tests";
    }

    virtual ReporterPreferences getPreferences() const {
      ReporterPreferences prefs;
      prefs.shouldRedirectStdOut = false;
      return prefs;
    }

    virtual void noMatchingTestCases(std::string const & spec) {
      if (OHMMS::Controller->comm.root()) {
        ConsoleReporter::noMatchingTestCases(spec);
       }
    }


    virtual void assertionStarting( AssertionInfo const &info) {
        ConsoleReporter::assertionStarting(info);
    }

    virtual bool assertionEnded(AssertionStats const& _assertionStats) {
      // prevent the intro message from being printed on non-root nodes
      if (!OHMMS::Controller->comm.root()) {
        currentTestRunInfo.used = true;
      }
      return ConsoleReporter::assertionEnded(_assertionStats);
    }

    virtual void testRunEnded(TestRunStats const& _testRunStats) {
      if (OHMMS::Controller->comm.root()) {
        ConsoleReporter::testRunEnded(_testRunStats);
      }
    }
  };
 
  CATCH_REGISTER_REPORTER("mpi", MPIReporter)

  MPIReporter::~MPIReporter() {}
}

#ifndef JetResolution_h
#define JetResolution_h

#include "JetResolutionObject.h"

#if 0
namespace edm {
    class EventSetup;
}
#endif


namespace JME {
    class JetResolution {
        public:
            JetResolution(const std::string& filename);
            JetResolution(const JetResolutionObject& object);
            JetResolution() {
                // Empty
            }

#if 0
            static const JetResolution get(const edm::EventSetup&, const std::string&);
#endif

            float getResolution(const JetParameters& parameters) const;

            void dump() const {
                m_object->dump();
            }

            // Advanced usage
            const JetResolutionObject* getResolutionObject() const {
                return m_object.get();
            }

        private:
            std::shared_ptr<JetResolutionObject> m_object;
    };

    class JetResolutionScaleFactor {
        public:
            JetResolutionScaleFactor(const std::string& filename);
            JetResolutionScaleFactor(const JetResolutionObject& object);
            JetResolutionScaleFactor() {
                // Empty
            }

#if 0
            static const JetResolutionScaleFactor get(const edm::EventSetup&, const std::string&);
#endif

            float getScaleFactor(const JetParameters& parameters, Variation variation = Variation::NOMINAL) const;

            void dump() const {
                m_object->dump();
            }

            // Advanced usage
            const JetResolutionObject* getResolutionObject() const {
                return m_object.get();
            }

        private:
            std::shared_ptr<JetResolutionObject> m_object;
    };

}

#endif


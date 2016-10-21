#ifndef SETUP_PARSER_H
#define SETUP_PARSER_H

#include "rapidxml.hpp"
#include "libfreenect2_data_structures.h"
#include <string>
#include <vector>

class SetupParser
{
    public:
        SetupParser();
        ~SetupParser();
        bool init(std::string xmlfilename);        
        std::vector<Parameters> getParameters();
        std::string getSetupName();
    private:
        void init_parameters();
        void insertParameterValue(Parameters &params, rapidxml::xml_node<>* node);

        char* text_;
        std::string setup_name_;
        rapidxml::xml_document<> doc_;
        std::vector<Parameters> parameter_vec_;

};


#endif


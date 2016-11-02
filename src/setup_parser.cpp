#include "setup_parser.h"
#include <fstream>
#include <cstdlib>
#include <iostream> 
#include <stdio.h>
#include <string.h>

SetupParser::SetupParser()
{

}

SetupParser::~SetupParser()
{
    delete text_;
}


bool SetupParser::init(std::string xmlfilename)
{
    std::ifstream xmlfile(xmlfilename.c_str(), std::ios_base::in);
    
    if(xmlfile)
    {
        xmlfile.seekg(0, xmlfile.end);
        int length = xmlfile.tellg();
        xmlfile.seekg(0, xmlfile.beg);
        
        text_ = new char[length];
        xmlfile.read(text_,length);
        doc_.parse<0>(text_);
        xmlfile.close();
        init_parameters();
        return true;
    }
    else
        return false;
}

void SetupParser::init_parameters()
{
    rapidxml::xml_node<>* root_node = doc_.first_node("setup");
    setup_name_ = root_node->first_attribute("name")->value();
    
    rapidxml::xml_node<>* current_node = root_node->first_node("pipeline");

    //loops over all pipelines
    while(current_node != 0)
    {
        Parameters tmp_param;
        tmp_param.setup_name = current_node->first_attribute("setup_name")->value();
        tmp_param.pipeline = current_node->first_attribute("name")->value();
        
        rapidxml::xml_node<>* current_parameter_node = current_node->first_node()->first_node();
        while(current_parameter_node != 0)
        {
            insertParameterValue(tmp_param, current_parameter_node);
            current_parameter_node = current_parameter_node->next_sibling();
        }
        
        parameter_vec_.push_back(tmp_param);
        current_node = current_node->next_sibling("pipeline");
    }

}

void SetupParser::insertParameterValue(Parameters &params, rapidxml::xml_node<>* node)
{
    //strcmp(str1,str2)
    if(strcmp(node->name(),"enable_bilateral_filter") == 0)
        params.enable_bilateral_filter = strcmp(node->value(), "1");
    else if(strcmp(node->name(),"enable_edge_filter") == 0)
        params.enable_edge_filter = strcmp(node->value(), "1");
    else if(strcmp(node->name(),"joint_bilateral_ab_threshold") == 0)
        params.joint_bilateral_ab_threshold = std::atof(node->value());
    else if(strcmp(node->name(),"joint_bilateral_max_edge") == 0)
        params.joint_bilateral_max_edge = std::atof(node->value());
    else if(strcmp(node->name(),"joint_bilateral_exp") == 0)
        params.joint_bilateral_exp = std::atof(node->value());
    else if(strcmp(node->name(),"gaussian_kernel[0]") == 0)
        params.gaussian_kernel[0] = std::atof(node->value());
    else if(strcmp(node->name(),"gaussian_kernel[1]") == 0)
        params.gaussian_kernel[1] = std::atof(node->value());
    else if(strcmp(node->name(),"gaussian_kernel[2]") == 0)
        params.gaussian_kernel[2] = std::atof(node->value());
    else if(strcmp(node->name(),"gaussian_kernel[3]") == 0)
        params.gaussian_kernel[3] = std::atof(node->value());
    else if(strcmp(node->name(),"gaussian_kernel[4]") == 0)
        params.gaussian_kernel[4] = std::atof(node->value());
    else if(strcmp(node->name(),"gaussian_kernel[5]") == 0)
        params.gaussian_kernel[5] = std::atof(node->value());
    else if(strcmp(node->name(),"gaussian_kernel[6]") == 0)
        params.gaussian_kernel[6] = std::atof(node->value());
    else if(strcmp(node->name(),"gaussian_kernel[7]") == 0)
        params.gaussian_kernel[7] = std::atof(node->value());
    else if(strcmp(node->name(),"gaussian_kernel[8]") == 0)
        params.gaussian_kernel[8] = std::atof(node->value());
    else if(strcmp(node->name(),"phase_offset") == 0)
        params.phase_offset = std::atof(node->value());
    else if(strcmp(node->name(),"unambigious_dist") == 0)
        params.unambigious_dist = std::atof(node->value());
    else if(strcmp(node->name(),"individual_ab_threshold") == 0)
        params.individual_ab_threshold = std::atof(node->value());
    else if(strcmp(node->name(),"ab_threshold") == 0)
        params.ab_threshold = std::atof(node->value());
    else if(strcmp(node->name(),"ab_confidence_slope") == 0)
        params.ab_threshold = std::atof(node->value());
    else if(strcmp(node->name(),"ab_confidence_offset") == 0)
        params.ab_confidence_offset = std::atof(node->value());
    else if(strcmp(node->name(),"min_dealias_confidence") == 0)
        params.min_dealias_confidence = std::atof(node->value());
    else if(strcmp(node->name(),"max_dealias_confidence") == 0)
        params.max_dealias_confidence = std::atof(node->value());
    else if(strcmp(node->name(),"edge_ab_avg_min_value") == 0)
        params.edge_ab_avg_min_value = std::atof(node->value());
    else if(strcmp(node->name(),"edge_ab_std_dev_threshold") == 0)
        params.edge_ab_std_dev_threshold = std::atof(node->value());
    else if(strcmp(node->name(),"edge_close_delta_threshold") == 0)
        params.edge_close_delta_threshold = std::atof(node->value());
    else if(strcmp(node->name(),"edge_far_delta_threshold") == 0)
        params.edge_far_delta_threshold = std::atof(node->value());
    else if(strcmp(node->name(),"edge_max_delta_threshold") == 0)
        params.edge_max_delta_threshold = std::atof(node->value());
    else if(strcmp(node->name(),"edge_avg_delta_threshold") == 0)
        params.edge_avg_delta_threshold = std::atof(node->value());
    else if(strcmp(node->name(),"max_edge_count") == 0)
        params.max_edge_count = std::atof(node->value());
    else if(strcmp(node->name(),"kde_sigma_sqr") == 0)
        params.kde_sigma_sqr = std::atof(node->value());
    else if(strcmp(node->name(),"unwrapping_likelihood_scale") == 0)
        params.unwrapping_likelihood_scale = std::atof(node->value());
    else if(strcmp(node->name(),"phase_confidence_scale") == 0)
        params.phase_confidence_scale = std::atof(node->value());
    else if(strcmp(node->name(),"kde_threshold") == 0)
        params.kde_threshold = std::atof(node->value());
    else if(strcmp(node->name(),"kde_neigborhood_size") == 0)
        params.kde_neigborhood_size = std::atoi(node->value());
    else if(strcmp(node->name(),"num_hyps") == 0)
        params.num_hyps = std::atoi(node->value());
    else if(strcmp(node->name(),"min_depth") == 0)
        params.min_depth = std::atof(node->value());
    else if(strcmp(node->name(),"max_depth") == 0)
        params.max_depth = std::atof(node->value());
    else if(strcmp(node->name(),"phase_noise_prediction_method") == 0)
        params.phase_noise_prediction_method = node->value();
		else 
        std::cout<<"parameter "<<node->name()<<" not defined or changeable\n";
     
   
}

std::vector<Parameters> SetupParser::getParameters()
{
    return parameter_vec_;
}

std::string SetupParser::getSetupName()
{
    return setup_name_;
}

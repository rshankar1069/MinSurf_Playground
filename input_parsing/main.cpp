#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <map>
#include <vector>
#include "atmsp.h"

int main(int argc, char* argv[])
{
    std::string filename,line;
    std::string Nx,Ny;
    std::string bottom,right,top,left;
    
    /* Defining a map container to store the constant names and their values */
    std::map<std::string,float> consts;
    std::map<std::string,float>::iterator constit;
    /* Defining a vector container to store the variable names */
    std::vector<std::string> vars;
    std::vector<std::string> tokens;

    if(argc<1)
    {
        std::cout << "Please pass a file to be read with the program !!" << std::endl;
        abort();
    }
    else
        filename = argv[1];

    std::ifstream inp_file(filename);
    if(inp_file.is_open())
    {
        while(std::getline(inp_file,line))
        {
            if (line[0] == '#' || line == "" || line[0] == ' ')
                continue;
            else
            {
                for(auto i = strtok(&line[0], " "); i != NULL; i = strtok(NULL, " "))
                    tokens.push_back(i);
            }
        }
        inp_file.close();
    }
    else
    {
        std::cout << "Unable to open file" << std::endl;
        abort();
    }
    
    for(int i=0; i<tokens.size(); i+=2)
    {
        if(tokens[i] == "Nx")
            Nx = std::stoi(tokens[i+1]);
        else if(tokens[i] == "Ny")
            Ny = std::stoi(tokens[i+1]);
        else if(tokens[i] == "BC.BOTTOM")
            bottom = tokens[i+1];
        else if(tokens[i] == "BC.RIGHT")
            right = tokens[i+1];
        else if(tokens[i] == "BC.TOP")
            top = tokens[i+1];
        else if(tokens[i] == "BC.LEFT")
            left = tokens[i+1];
    }
    
    /* Variable name storage */
    auto iterator = std::find(tokens.begin(),tokens.end(),"NUM_VARS");
    if (iterator != tokens.cend())
    {
        int pos = std::distance(tokens.begin(),iterator);
        int numvars = std::stoi(tokens[pos+1]);
        pos = pos+2;
        for(int i=0; i<2*numvars; i+=2)
            vars.push_back(tokens[pos+i+1]); 
    }
    else
    {
        std::cout << "Please enter some variables to be parsed by the expression parser !!" << std::endl;
        abort();
    }

    /* Constants name storage */
    iterator = std::find(tokens.begin(),tokens.end(),"NUM_CONSTS");
    if (iterator != tokens.cend())
    {
        int pos = std::distance(tokens.begin(),iterator);
        int numconsts = std::stoi(tokens[pos+1]);
        pos = pos+2;
        for(int i=0; i<2*numconsts; i+=2)
        {
            consts.insert(std::pair<std::string,float>(tokens[pos+i],std::stof(tokens[pos+i+1])));
            // constsNames.push_back(tokens[pos+i]);
            // constsVals.push_back(std::stof(tokens[pos+i+1]));
        } 
    }
    else
    {
        std::cout << "No external constants defined !!" << std::endl;
    }

    /* Using the equation parser */

    ATMSP<double> parser;
    ATMSB<double> byteCode;
    // ATMSP<double> parserBottom,parserRight,parserTop,parserLeft;
    // ATMSB<double> byteCodeBottom,byteCodeRight,byteCodeTop,byteCodeLeft;
    std::string varnames;

    std::vector<double> xval;

    int N = 100;
    double dt = 1.0/N;
    int i=0;

    while(i*dt <= N*dt)
    {
        xval.push_back(dt*i);
        i++;
    }

    for(constit=consts.begin(); constit!=consts.end(); constit++)
    {
        parser.addConstant(constit->first,constit->second);
    }

    for(int i=0; i<vars.size(); i++)
    {
        if(i!=vars.size()-1)
            varnames += vars[i]+",";
        else
            varnames += vars[i];
    }

    parser.parse(byteCode,bottom,varnames);
    i = 0;
    std::string writefile = "output.txt";
    std::ofstream out_file(writefile);
    while(i*dt <= 1.0)
    {
        byteCode.var[0] = xval[i];
        out_file << xval[i] << "\t" << byteCode.run() << std::endl;
        i++;
    }
    return 0;
}
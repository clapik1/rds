#ifndef RDS_ARGSPARSER_H
#define RDS_ARGSPARSER_H


#include <vector>
#include <string>

class argsParser {
public:
    argsParser (int &argc, char **argv);
    std::string getOption(const std::string &option) const;
    bool optionExists(const std::string &option) const;
private:
    std::vector<std::string> tokens;
};


#endif //RDS_ARGSPARSER_H

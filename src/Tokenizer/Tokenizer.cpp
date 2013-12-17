/*
 *  This file is part of AQUAgpusph, a free CFD program based on SPH.
 *  Copyright (C) 2012  Jose Luis Cercos Pita <jl.cercos@upm.es>
 *
 *  AQUAgpusph is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  AQUAgpusph is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with AQUAgpusph.  If not, see <http://www.gnu.org/licenses/>.
 */

// ----------------------------------------------------------------------------
// Include the main header
// ----------------------------------------------------------------------------
#include <Tokenizer/Tokenizer.h>

Tokenizer::Tokenizer()
{
	// Register default variables
	registerVariable("pi", M_PI);
	registerVariable("e", M_E);
}

Tokenizer::~Tokenizer()
{
	unsigned int i;
	for(i=0;i<regNames.size();i++){
	    delete[] regNames.at(i);
	}
	regNames.clear();
	regValues.clear();
}

bool Tokenizer::registerVariable(const char* name, float value)
{
	int i = findVariable(name);
	if(i == -1){
	    char *Name = new char[32];
	    strcpy(Name, name);
	    regNames.push_back(Name);
	    regValues.push_back(value);
	    return false;
	}
	else{
	    regValues.at(i) = value;
	    return true;
	}
}

bool Tokenizer::unregisterVariable(const char* name)
{
	int i = findVariable(name);
	if(i == -1){
	    return false;
	}
	else{
	    delete[] regNames.at(i);
	    regNames.erase(regNames.begin() + i);
	    regValues.erase(regValues.begin() + i);
	    return true;
	}
}

float Tokenizer::variable(const char* name, bool whole, unsigned int *readed)
{
	int i = findVariable(name, whole, readed);
	if(i==-1)
	    return 0.f;
	return regValues.at(i);
}

int Tokenizer::findVariable(const char* name, bool whole, unsigned int *readed)
{
	unsigned int i;
	size_t n = strlen(name);
	for(i=0;i<regNames.size();i++){
	    if(!whole)
	        n = strlen(regNames.at(i));
	    if(!strncmp(name,regNames.at(i),n)){
	        if(readed)
	            *readed = n;
	        return (int)i;
	    }
	}
	if(readed){
	    *readed=0;
	}
	return -1;
}

bool Tokenizer::isVariable(const char* name, bool whole, unsigned int *readed)
{
	unsigned int i;
	size_t n = strlen(name);
	for(i=0;i<regNames.size();i++){
	    if(!whole)
	        n = strlen(regNames.at(i));
	    if(!strncmp(name,regNames.at(i),n)){
	        if(readed)
	            *readed = n;
	        return true;
	    }
	}
	if(readed){
	    *readed=0;
	}
	return false;
}

bool Tokenizer::isFloat(const char* word, bool whole, unsigned int *readed)
{
	if(readed)
	    *readed=0;
	unsigned int i,i0=0;
	char c = word[i0];
	if(c=='-'){
	    i0++;
	    c = word[i0];
	}
	if(! ( (c=='0') || (c=='1') || (c=='2') || (c=='3') || (c=='4') || (c=='5') || (c=='6') || (c=='7') || (c=='8') || (c=='9') ) ){
	    return false;
	}
	bool decimals=false;
	for(i=i0+1;i<strlen(word);i++){
	    c = word[i];
	    if(c=='.'){
	        if(decimals){
	            if(!whole){
	                if(readed)
	                    *readed = i;
	                return true;
	            }
	            else
	                return false;
	        }
	        decimals = true;
	        continue;
	    }
	    if(! ( (c=='0') || (c=='1') || (c=='2') || (c=='3') || (c=='4') || (c=='5') || (c=='6') || (c=='7') || (c=='8') || (c=='9') ) ){
	        if(!whole){
	            if(readed)
	                *readed = i;
	            return true;
	        }
	        else
	            return false;
	    }
	}
	// Whole word is a number
	if(readed)
	    *readed = strlen(word);
	return true;
}

bool Tokenizer::solve(const char* eq, float *value)
{
	// Re-format sentence
	char sentence[256];
	strcpy(sentence, format(eq));
	// Solve
	return _solve(sentence,value);
}

bool Tokenizer::_solve(const char* eq, float *value)
{
	// Test if is a whole number
	unsigned int readed=0;
	float aux;
	if(isFloat(eq,false,&readed)){
	    aux = atof(eq);
	    if(readed < strlen(eq)){                    // Remains equation to solve
	        if(operate(aux,&eq[readed],value)){    // That must continue with an operator
	            return true;
	        }
	        else{
	            return false;
	        }
	    }
	    *value = aux;
	    return true;
	}
	// Test if is a variable
	if(isVariable(eq,false,&readed)){
	    aux = variable(eq,false);
	    if(readed < strlen(eq)){                    // Remains equation to solve
	        if(operate(aux,&eq[readed],value)){      // That must continue with an operator
	            return true;
	        }
	        else{
	            return false;
	        }
	    }
	    *value = aux;
	    return true;
	}
	// Test if is a function
	if(solveFunction(eq, &aux, &readed)){
	    if(readed < strlen(eq)){                    // Remains equation to solve
	        if(operate(aux,&eq[readed],value)){    // That must continue with an operator
	            return true;
	        }
	        else{
	            return false;
	        }
	    }
	    *value = aux;
	    return true;
	}
	printf("\tERROR (Tokenizer::solve): Unhandled expression.\n");
	printf("\t\t\"%s\"\n", eq);
	return false;
}

const char* Tokenizer::format(const char* eq)
{
	static char arg[256];
	// Remove the spaces
	const char* aux = removeSpaces(eq);
	if(!aux)
	    return 0;
	strcpy(arg,aux);
	// Insert Parentheses
	aux = sortSentence(arg);
	if(!aux)
	    return 0;
	strcpy(arg,aux);
	return arg;
}

const char* Tokenizer::removeSpaces(const char* eq)
{
	static char arg[256]; strcpy(arg, eq);
	char *ReplacePoint=0;
	ReplacePoint = strstr(arg," ");
	while(ReplacePoint){
	    char StrBackup[256];
	    strcpy(StrBackup, &(ReplacePoint[1]));
	    strcpy(&(ReplacePoint[0]),StrBackup);
	    ReplacePoint = strstr(arg," ");
	}
	// Deletes also the line breaks
	ReplacePoint = strstr(arg,"\n");
	while(ReplacePoint){
	    char StrBackup[256];
	    strcpy(StrBackup, &(ReplacePoint[1]));
	    strcpy(&(ReplacePoint[0]),StrBackup);
	    ReplacePoint = strstr(arg,"\n");
	}
	return arg;
}

const char* Tokenizer::sortSentence(const char* eq)
{
	static char arg[256]; strcpy(arg, eq);
	char *oper=0;
	const char *aux=0;
	int i,count, n;
	// Look for operators decreasing priority
	oper = strstr(arg,"^");
	count = 0;
	while(oper){
	    count++;
	    n = strlen(arg) - strlen(oper);
	    aux = joinTerms(arg, n);
	    if(!aux)
	        return 0;
	    strcpy(arg, aux);
	    oper = strstr(arg,"^");
	    for(i=0;i<count;i++)
	        oper = strstr(&oper[1],"^");
	}
	oper = strstr(arg,"*");
	count = 0;
	while(oper){
	    count++;
	    n = strlen(arg) - strlen(oper);
	    aux = joinTerms(arg, n);
	    if(!aux)
	        return 0;
	    strcpy(arg, aux);
	    oper = strstr(arg,"*");
	    for(i=0;i<count;i++)
	        oper = strstr(&oper[1],"*");
	}
	oper = strstr(arg,"/");
	count = 0;
	while(oper){
	    count++;
	    n = strlen(arg) - strlen(oper);
	    aux = joinTerms(arg, n);
	    if(!aux)
	        return 0;
	    strcpy(arg, aux);
	    oper = strstr(arg,"/");
	    for(i=0;i<count;i++)
	        oper = strstr(&oper[1],"/");
	}
	oper = strstr(arg,"+");
	count = 0;
	while(oper){
	    count++;
	    n = strlen(arg) - strlen(oper);
	    aux = joinTerms(arg, n);
	    if(!aux)
	        return 0;
	    strcpy(arg, aux);
	    oper = strstr(arg,"+");
	    for(i=0;i<count;i++)
	        oper = strstr(&oper[1],"+");
	}
	oper = strstr(arg,"-");
	count = 0;
	while(oper){
	    count++;
	    n = strlen(arg) - strlen(oper);
	    aux = joinTerms(arg, n);
	    if(!aux)
	        return 0;
	    strcpy(arg, aux);
	    oper = strstr(arg,"-");
	    for(i=0;i<count;i++)
	        oper = strstr(&oper[1],"-");
	}
	return arg;
}

const char* Tokenizer::joinTerms(const char* eq, int n)
{
	static char arg[256]; strcpy(arg, eq);
	int i;
	int balance=0, lastPar=0;
	// Look for start of group
	for(i=n-1;i>=0;i--){
	    // Parentheses
	    if(arg[i] == ')'){
	        balance++;
	        if(i>lastPar){
	            lastPar = i;
	        }
	    }
	    else if(arg[i] == '('){
	        balance--;
	        if(balance <= 0){
	            strcpy(&arg[i+1],"(");
	            strcpy(&arg[i+2],&eq[i+1]);
	            break;
	        }
	    }
	    // Change of argument
	    if(arg[i] == ','){
	        if(balance <= 0){
	            strcpy(&arg[i+1],"(");
	            strcpy(&arg[i+2],&eq[i+1]);
	            break;
	        }
	    }
	    // Operators
	    if(balance <= 0){
	        // Special negative number
	        if( arg[i] == '-' ){
	            if(i==0){       // Negative number at the start of the sentence (We use as start of sentence reached)
	                i=-1;
	                break;
	            }
	            else if((arg[i-1] == '(') || (arg[i-1] == ',')){    // Negative number at the start of argument
	                strcpy(&arg[i],"(");
	                strcpy(&arg[i+1],&eq[i]);
	                break;
	            }
	        }
	        if( (arg[i] == '+') || (arg[i] == '-') || (arg[i] == '*') || (arg[i] == '/') || (arg[i] == '^') ){
	            strcpy(&arg[i+1],"(");
	            strcpy(&arg[i+2],&eq[i+1]);
	            break;
	        }
	    }
	}
	if(balance > 0){
	    printf("\tERROR (Tokenizer::joinTerms): Unbalanced parentheses.\n");
	    printf("\t\t\"%s\"\n\t\t", arg);
	    for(i=0;i<lastPar;i++){
	        printf(" ");
	    }
	    printf("^\n");
	    return 0;
	}
	if(i<0){        // Start of sentence reached
	        strcpy(&arg[0],"(");
	        strcpy(&arg[1],&eq[0]);
	}
	// Look for end of group
	balance = 0;
	for(i=n+2;i<(int)strlen(arg);i++){       // n+2 because we need jump over the operator, and count the inserted parentheses
	    // Parentheses
	    if(arg[i] == '('){
	        balance++;
	        if(i>lastPar)
	            lastPar = i;
	    }
	    else if(arg[i] == ')'){
	        balance--;
	        if(balance <= 0){
	            strcpy(&arg[i],")");
	            strcpy(&arg[i+1],&eq[i-1]);
	            break;
	        }
	    }
	    // Change of argument
	    if(arg[i] == ','){
	        if(balance <= 0){
	            strcpy(&arg[i],")");
	            strcpy(&arg[i+1],&eq[i-1]);
	            break;
	        }
	    }
	    // Operators
	    if(balance <= 0){
	        if( (arg[i] == '+') || (arg[i] == '-') || (arg[i] == '*') || (arg[i] == '/') || (arg[i] == '^') ){
	            strcpy(&arg[i],")");
	            strcpy(&arg[i+1],&eq[i-1]);
	            break;
	        }
	    }
	}
	if(balance > 0){
	    printf("\tERROR (Tokenizer::joinTerms): Unbalanced parentheses.\n");
	    printf("\t\t\"%s\"\n\t\t", arg);
	    for(i=0;i<lastPar;i++)
	        printf(" ");
	    printf("^\n");
	    return 0;
	}
	if(i==(int)strlen(arg)){        // End of sentence reached
	        strcpy(&arg[i],")");
	}
	return arg;
}

bool Tokenizer::operate(float leftOp, const char* op, float *value)
{
	//! Solve right operator
	float rightOp = 0.0;
	if(!_solve(&op[1],&rightOp))
	    return false;
	//! Solve left vs. right operation
	if(op[0] == '+'){    // Add
	    *value = leftOp + rightOp;
	    return true;
	}
	if(op[0] == '-'){    // Substract
	    *value = leftOp - rightOp;
	    return true;
	}
	if(op[0] == '*'){    // Multiply
	    *value = leftOp * rightOp;
	    return true;
	}
	if(op[0] == '/'){    // Divide
	    *value = leftOp / rightOp;
	    return true;
	}
	if(op[0] == '^'){    // Pow
	    *value = pow(leftOp,rightOp);
	    return true;
	}
	printf("\tERROR (Tokenizer::operate): Any valid operator found.\n");
	printf("\t\t\"%s\"\n", op);
	return false;
}

bool Tokenizer::solveFunction(const char* func, float *value, unsigned int *readed)
{
	//! Look for know functions
	// Parentheses "function"
	char eq[256];
	if(!strncmp(func, "(", 1)){
	    strcpy(eq,func);
	    if(solveParentheses(eq,value,readed)){
	        return true;
	    }
	    return false;
	}
	else if(!strncmp(func, "-", 1)){
	    strcpy(eq,&func[1]);
	    if(solveNegative(eq,value,readed)){
	        if(readed)
	            *readed += 1;   // sin word is not taken account
	        return true;
	    }
	    return false;
	}
	else if(!strncmp(func, "rad", 3)){
	    strcpy(eq,&func[3]);
	    if(solveRad(eq,value,readed)){
	        if(readed)
	            *readed += 3;   // sin word is not taken account
	        return true;
	    }
	    return false;
	}
	else if(!strncmp(func, "deg", 3)){
	    strcpy(eq,&func[3]);
	    if(solveDeg(eq,value,readed)){
	        if(readed)
	            *readed += 3;   // sin word is not taken account
	        return true;
	    }
	    return false;
	}
	else if(!strncmp(func, "sin", 3)){
	    strcpy(eq,&func[3]);
	    if(solveSin(eq,value,readed)){
	        if(readed)
	            *readed += 3;   // sin word is not taken account
	        return true;
	    }
	    return false;
	}
	else if(!strncmp(func, "cos", 3)){
	    strcpy(eq,&func[3]);
	    if(solveCos(eq,value,readed)){
	        if(readed)
	            *readed += 3;   // sin word is not taken account
	        return true;
	    }
	    return false;
	}
	else if(!strncmp(func, "tan", 3)){
	    strcpy(eq,&func[3]);
	    if(solveTan(eq,value,readed)){
	        if(readed)
	            *readed += 3;   // sin word is not taken account
	        return true;
	    }
	    return false;
	}
	else if(!strncmp(func, "asin", 4)){
	    strcpy(eq,&func[4]);
	    if(solveAsin(eq,value,readed)){
	        if(readed)
	            *readed += 4;   // sin word is not taken account
	        return true;
	    }
	    return false;
	}
	else if(!strncmp(func, "acos", 4)){
	    strcpy(eq,&func[4]);
	    if(solveAcos(eq,value,readed)){
	        if(readed)
	            *readed += 4;   // sin word is not taken account
	        return true;
	    }
	    return false;
	}
	else if(!strncmp(func, "atan", 4)){
	    strcpy(eq,&func[4]);
	    if(solveAtan(eq,value,readed)){
	        if(readed)
	            *readed += 4;   // sin word is not taken account
	        return true;
	    }
	    return false;
	}
	else if(!strncmp(func, "atan2", 5)){
	    strcpy(eq,&func[5]);
	    if(solveAtan2(eq,value,readed)){
	        if(readed)
	            *readed += 5;   // sin word is not taken account
	        return true;
	    }
	    return false;
	}
	else if(!strncmp(func, "exp", 3)){
	    strcpy(eq,&func[3]);
	    if(solveExp(eq,value,readed)){
	        if(readed)
	            *readed += 3;   // sin word is not taken account
	        return true;
	    }
	    return false;
	}
	else if(!strncmp(func, "log", 3)){
	    strcpy(eq,&func[3]);
	    if(solveLog(eq,value,readed)){
	        if(readed)
	            *readed += 3;   // sin word is not taken account
	        return true;
	    }
	    return false;
	}
	else if(!strncmp(func, "log10", 5)){
	    strcpy(eq,&func[5]);
	    if(solveLog(eq,value,readed)){
	        if(readed)
	            *readed += 5;   // sin word is not taken account
	        return true;
	    }
	    return false;
	}
	else if(!strncmp(func, "pow", 3)){
	    strcpy(eq,&func[3]);
	    if(solvePow(eq,value,readed)){
	        if(readed)
	            *readed += 3;   // sin word is not taken account
	        return true;
	    }
	    return false;
	}
	else if(!strncmp(func, "sqrt", 4)){
	    strcpy(eq,&func[4]);
	    if(solveSqrt(eq,value,readed)){
	        if(readed)
	            *readed += 4;   // sin word is not taken account
	        return true;
	    }
	    return false;
	}
	else if(!strncmp(func, "clamp", 5)){
	    strcpy(eq,&func[5]);
	    if(solveClamp(eq,value,readed)){
	        if(readed)
	            *readed += 5;   // sin word is not taken account
	        return true;
	    }
	    return false;
	}
	else if(!strncmp(func, "min", 3)){
	    strcpy(eq,&func[3]);
	    if(solveMin(eq,value,readed)){
	        if(readed)
	            *readed += 3;   // sin word is not taken account
	        return true;
	    }
	    return false;
	}
	else if(!strncmp(func, "max", 3)){
	    strcpy(eq,&func[3]);
	    if(solveMax(eq,value,readed)){
	        if(readed)
	            *readed += 3;   // sin word is not taken account
	        return true;
	    }
	    return false;
	}
	return false;
}

unsigned int Tokenizer::getNumberOfArguments(const char* func)
{
	const char *separator=0;
	separator = strstr(func,",");
	unsigned int i, count = 1;
	while(separator){
	    count++;
	    separator = strstr(func,",");
	    for(i=1;i<count;i++)
	        separator = strstr(&separator[1],",");
	}
	return count;
}

const char* Tokenizer::getArgument(const char* func, unsigned int id)
{
	static char arg[256]; strcpy(arg,"");
	const char *separator=0;
	unsigned int i, count = 0, n;
	int balance;
	if(id==0){
	    strcpy(arg,&func[1]);
	}
	else{
	    separator = func;
	    while(separator){
	        separator = strstr(&separator[1],",");
	        if(!separator)
	            break;
	        n = strlen(func) - strlen(separator);
	        // Look for balanced parentheses (if unbalanced parentheses, can be into a function)
	        balance = 0;
	        for(i=0;i<n;i++){
	            if(func[i] == '(')
	               balance++;
	            else if(func[i] == ')')
	                balance--;
	        }
	        if(balance != 1){       // Open parenthes must be taken into account
	            continue;
	        }
	        count++;
	        if(count==id){
	            strcpy(arg,&separator[1]);
	            break;
	        }
	    }
	}

	if(count < id){
	    printf("\tERROR (Tokenizer::getArgument): Not enought arguments.\n");
	    printf("\t\t%u request, %u encountered.\n", id, count);
	    return 0;
	}
	// Look for the end of the argument
	separator = strstr(arg,",");
	if(separator){
	    while(separator){
	        n = strlen(arg) - strlen(separator);
	        // Look for balanced parentheses (if unbalanced parentheses, can be into a function)
	        balance = 0;
	        for(i=0;i<n;i++){
	            if(arg[i] == '(')
	               balance++;
	            else if(arg[i] == ')')
	                balance--;
	        }
	        if(balance == 0){
	            break;
	        }
	        separator = strstr(&separator[1],",");
	    }
	    if(balance != 0){ // Any other separator found, must be the end of the expression
	        n = strlen(arg) - 1;
	    }
	}
	else
	    n = strlen(arg) - 1;
	strcpy(&arg[n], "");
	return arg;
}

bool Tokenizer::solveParentheses(const char* func, float *value, unsigned int *readed){
	char arg[256]; strcpy(arg, getParenthesesArg(func, readed));
	if(!_solve(arg,value)){
	    return false;
	}
	return true;
}

const char* Tokenizer::getParenthesesArg(const char* eq, unsigned int *readed)
{
	if(readed)
	    *readed = 0;
	if(eq[0] != '('){
	    printf("\tERROR (Tokenizer::getParenthesesArg): Expression don't start with parentheses.\n");
	    printf("\t\t\"%s\"\n", eq);
	    return 0;
	}
	unsigned int i;
	unsigned int balance=1;
	for(i=1;i<strlen(eq);i++){
	    if(eq[i] == '(')
	        balance++;
	    else if(eq[i] == ')')
	        balance--;
	    if(!balance){           // Parentheses are balanced, we get the string selected
	        if(i==1){
	            printf("\tERROR (Tokenizer::getParenthesesArg): Empty expression.\n");
	            printf("\t\t\"%s\"\n", eq);
	            return 0;
	        }
	        static char arg[256];
	        strncpy(arg,&eq[1],i-1);strcpy(&arg[i-1], "");
	        if(readed)
	            *readed = i+1;
	        return arg;
	    }
	}
	printf("\tERROR (Tokenizer::getParenthesesArg): Parentheses unbalanced.\n");
	printf("\t\t\"%s\"\n",eq);
	return 0;
}

bool Tokenizer::solveNegative(const char* func, float *value, unsigned int *readed){
	char arg[256]; strcpy(arg, getParenthesesArg(func, readed));

	if(!_solve(arg,value))
	    return false;
	*value = -*value;
	return true;
}

bool Tokenizer::solveRad(const char* func, float *value, unsigned int *readed){
	char arg[256]; strcpy(arg, getParenthesesArg(func, readed));

	if(!_solve(arg,value))
	    return false;
	*value = *value * M_PI / 180.f;
	return true;
}

bool Tokenizer::solveDeg(const char* func, float *value, unsigned int *readed){
	char arg[256]; strcpy(arg, getParenthesesArg(func, readed));

	if(!_solve(arg,value))
	    return false;
	*value = *value * 180.f / M_PI;
	return true;
}

bool Tokenizer::solveSin(const char* func, float *value, unsigned int *readed){
	char arg[256]; strcpy(arg, getParenthesesArg(func, readed));

	if(!_solve(arg,value))
	    return false;
	*value = sin(*value);
	return true;
}

bool Tokenizer::solveCos(const char* func, float *value, unsigned int *readed){
	char arg[256]; strcpy(arg, getParenthesesArg(func, readed));

	if(!_solve(arg,value))
	    return false;
	*value = cos(*value);
	return true;
}

bool Tokenizer::solveTan(const char* func, float *value, unsigned int *readed){
	char arg[256]; strcpy(arg, getParenthesesArg(func, readed));

	if(!_solve(arg,value))
	    return false;
	*value = tan(*value);
	return true;
}

bool Tokenizer::solveAsin(const char* func, float *value, unsigned int *readed){
	char arg[256]; strcpy(arg, getParenthesesArg(func, readed));

	if(!_solve(arg,value))
	    return false;
	*value = asinf(*value);
	return true;
}

bool Tokenizer::solveAcos(const char* func, float *value, unsigned int *readed){
	char arg[256]; strcpy(arg, getParenthesesArg(func, readed));

	if(!_solve(arg,value))
	    return false;
	*value = acosf(*value);
	return true;
}

bool Tokenizer::solveAtan(const char* func, float *value, unsigned int *readed){
	char arg[256]; strcpy(arg, getParenthesesArg(func, readed));

	if(!_solve(arg,value))
	    return false;
	*value = atanf(*value);
	return true;
}

bool Tokenizer::solveAtan2(const char* func, float *value, unsigned int *readed){
	char arg[256]; strcpy(arg, getParenthesesArg(func, readed));    // Only to get readed, and read for errors
	//! Get y argument
	float arg1;
	strcpy(arg, "");
	strcpy(arg, getArgument(func, 0));
	if(!_solve(arg,&arg1))
	    return false;
	//! Get x argument
	float arg2;
	strcpy(arg, "");
	strcpy(arg, getArgument(func, 1));
	if(!_solve(arg,&arg2))
	    return false;
	//! Calculate solution
	*value = atan2(arg1,arg2);
	return true;
}

bool Tokenizer::solveExp(const char* func, float *value, unsigned int *readed){
	char arg[256]; strcpy(arg, getParenthesesArg(func, readed));

	if(!_solve(arg,value))
	    return false;
	*value = exp(*value);
	return true;
}

bool Tokenizer::solveLog(const char* func, float *value, unsigned int *readed){
	char arg[256]; strcpy(arg, getParenthesesArg(func, readed));

	if(!_solve(arg,value))
	    return false;
	*value = log(*value);
	return true;
}

bool Tokenizer::solveLog10(const char* func, float *value, unsigned int *readed){
	char arg[256]; strcpy(arg, getParenthesesArg(func, readed));

	if(!_solve(arg,value))
	    return false;
	*value = log10(*value);
	return true;
}

bool Tokenizer::solvePow(const char* func, float *value, unsigned int *readed){
	char arg[256]; strcpy(arg, getParenthesesArg(func, readed));    // Only to get readed, and read for errors
	//! Get argument 1
	float arg1;
	strcpy(arg, "");
	strcpy(arg, getArgument(func, 0));
	if(!_solve(arg,&arg1))
	    return false;
	//! Get argument 2
	float arg2;
	strcpy(arg, "");
	strcpy(arg, getArgument(func, 1));
	if(!_solve(arg,&arg2))
	    return false;
	//! Calculate solution
	*value = pow(arg1,arg2);
	return true;
}

bool Tokenizer::solveSqrt(const char* func, float *value, unsigned int *readed){
	char arg[256]; strcpy(arg, getParenthesesArg(func, readed));

	if(!_solve(arg,value))
	    return false;
	*value = sqrt(*value);
	return true;
}

bool Tokenizer::solveClamp(const char* func, float *value, unsigned int *readed){
	char arg[256]; strcpy(arg, getParenthesesArg(func, readed));    // Only to get readed, and read for errors
	//! Get argument 1 (value to adjust)
	float arg1;
	strcpy(arg, getArgument(func, 0));
	if(!_solve(arg,&arg1))
	    return false;
	//! Get argument 2 (minimum value)
	float arg2;
	strcpy(arg, getArgument(func, 1));
	if(!_solve(arg,&arg2))
	    return false;
	//! Get argument 3 (maximum value)
	float arg3;
	strcpy(arg, getArgument(func, 2));
	if(!_solve(arg,&arg3))
	    return false;
	//! Calculate solution
	*value = clamp(arg1,arg2,arg3);
	return true;
}

bool Tokenizer::solveMin(const char* func, float *value, unsigned int *readed){
	char arg[256]; strcpy(arg, getParenthesesArg(func, readed));    // Only to get readed, and read for errors
	//! Get argument 1
	float arg1;
	strcpy(arg, getArgument(func, 0));
	if(!_solve(arg,&arg1))
	    return false;
	//! Get argument 2
	float arg2;
	strcpy(arg, getArgument(func, 1));
	if(!_solve(arg,&arg2))
	    return false;
	//! Calculate solution
	*value = min(arg1,arg2);
	return true;
}

bool Tokenizer::solveMax(const char* func, float *value, unsigned int *readed){
	char arg[256]; strcpy(arg, getParenthesesArg(func, readed));    // Only to get readed, and read for errors
	//! Get argument 1
	float arg1;
	strcpy(arg, getArgument(func, 0));
	if(!_solve(arg,&arg1))
	    return false;
	//! Get argument 2
	float arg2;
	strcpy(arg, getArgument(func, 1));
	if(!_solve(arg,&arg2))
	    return false;
	//! Calculate solution
	*value = max(arg1,arg2);
	return true;
}

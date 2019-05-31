#ifndef LIB_H
#define LIB_H


// headers
#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <strstream>
#include <cctype>

// definitions
#define MLL_lib 1048576

// structures
//using std::ifstream;
using namespace std;
struct charplus
{
	char*** x;
	int* m;
	int n;
};

struct intplus
{
	int** x;
	int* m;
	int n;
};

struct floatplus
{
//	float** x;
	// ---data from body of file.
	double** x;
	int* m;
	int n;

	// ---data from comments.
	int ndat;
	double* ddat;
};

struct plexon
{
	double* t;
	int* c;
	int n;
};

// prototypes

char* cat2(char* prefix, char* suffix);

int choose(int k, int n);

double drand();
double drand(float x);
double drand(float xmin, float xmax);

double entropy(int n, int* p);
double entropy(int n, float* p, int norm=0);
double entropy(int n, double* p, int norm=0);

double entropy2(int n, int* p);
double entropy2(int n, float* p, int norm);
double entropy2(int n, double* p, int norm);

charplus get_chars(FILE* f);
float** get_cols(int n_col, int* col, int& n, FILE* f, int silent=0);

int** get_cols(
	FILE* f, int& lines, int& n_elements, int flag, int min_elements=0);
float** get_cols(
	FILE* f, int& lines, int& n_elements, float flag, int min_elements=0);
double** get_cols(
	FILE* f, int& lines, int& n_elements, double flag, int min_elements=0);

floatplus get_cols(FILE* f);
floatplus get_cols(
	FILE* f, int ndat, char** dat, double dflt, char sep,
	int min_elements=1);

int* header_info(FILE* f, int n, char** s, int dflt);

double* parse_file(
	FILE* f, ostrstream& x_ptr, ostrstream& m_ptr, int& lines,
	int ndat, char** dat, double dflt, char sep, int min_elements);
int xdat(int ndat, int* cdat, char** dat, char sep,
	int n_line, int* nw, char* line, double* ddat);
int isnumber(char* s, double& x);
double kl_divergence(int n, float* p, float* q);
int length_string(char* s);

int max(int n, int* x);
float max(int n, float* x);
double max(int n, double* x);
int min(int n, int* x);
float min(int n, float* x);
double min(int n, double* x);

double* mmult(int n, int m, double** a, double* x);
double** mmult(int n, int m, int l, double** a, double** x);


char** cchar(int n1, int n2);
int** cint(int n1, int n2);
float** cfloat(int n1, int n2);
double** cdouble(int n1, int n2);
double** cdouble(int n1, int* n2);
void cfree(int n, char** x);
void cfree(int n, int** x);
void cfree(int n, float** x);
void cfree(int n, double** x);

char** newchar(int n1, int n2);
double** newdouble(int n1, int n2);
double** newdouble(int n1, int* n2);

float** newfloat(int n1, int n2);
float** newfloat(int n1, int* n2);

int** newint(int n1, int n2);
int** newint(int n1, int* n2);

int* init_constant(int n, int constant);
int** init_constant(int n, int m, int constant);
int** init_constant(int n, int* m, int constant);
float* init_constant(int n, float constant);
float** init_constant(int n, int m, float constant);
float** init_constant(int n, int* m, float constant);
double* init_constant(int n, double constant);
double** init_constant(int n, int m, double constant);
double** init_constant(int n, int* m, double constant);


float mi(float** p, int imax, int jmax);

int rabs(int x);
float rabs(float x);
double rabs(double x);

unsigned short* read_binary(FILE* f, int& n, short dum);
int* read_binary(FILE* f, int& n, int dum);
float* read_binary(FILE* f, int& n, float dum);
double* read_binary(FILE* f, int& n, double dum);
int* read_binary(char* infile, int& n, int dum);
float* read_binary(char* infile, int& n, float dum);
double* read_binary(char* infile, int& n, double dum);

double* read_n(FILE* f, int lines, int col, int& n, int silent);
plexon read_plexon(FILE* f, FILE* g, int sort);

int readline(ifstream& infile, int nmax, int line, int& nchars, char* buf);
int* parse_c(char* buf, int& n);

int parse_comline(int n, char** t, char* s, int inc, int dflt);
float parse_comline(int n, char** t, char* s, int inc, float dflt);
double parse_comline(int n, char** t, char* s, int inc, double dflt);
int* parse_comline(int n, char** t, char* s, int& nmatch);
int* parse_comline(int n, char** t, char* s, int dum, int& nmatch);
int* parse_comline(int n, char** t, char* s, int dum, int& nmatch, char* dflt);
float* parse_comline(
	int n, char** t, char* s, float dum, int& nmatch, char* dflt);
double* parse_comline(
	int n, char** t, char* s, double dum, int& nmatch, char* dflt);
char* parse_comline(int n, char** t, char* s, int inc, char* dflt);
double* s_parse(int n, char* v);

int unit_no(char* s);

void write_err(char* s);
void write_err(char* s, int n);
void write_err(char* s1, int n, char* s2);
void write_err(char* s, double z);
void write_err(char* s1, char* s2);
void write_err(char* s1, char* s2, char* s3);
void write_err(char* s1, char* s2, int n, char* s3);
void write_err(char* s1, char* s2, int n1, char* s3, int n2);
void write_err(char* s1, char* s2, char* s3, int n, char* s4);
void write_err(char* s1, char* s2, char* s3, int n1, char* s4, int n2);
void write_err(char* s1, char* s2, char* s3, char* s4, char* s5);
void write_warning(char* s);
void write_warning(char* s1, char* s2, char* s3);
void write_warning(char* s1, int n, char* s2);
void write_warning(char* s1, int n1, char* s2, int n2, char* s3);
void write_warning(char* s1, char* s2, char* s3, int n, char* s4);
void write_warning(char* s1, int n1, char* s2, char* s3, int n2, char* s4);
void write_warning(char* s1, double x, char* s2);
void write_message(char* s, char* t);

void zfree(floatplus& z);

#endif

/* omdlib
 ********************************************************
 *name
 * ObjectMD
 * Molecular Dynamics Class Library
 * Version 2.0 (2009)
 *
 * Released under:
 *       GNU GENERAL PUBLIC LICENSE Version 3
 *       Copyright (C) 2005  Yudi Rosandi
 *
 * Read license.txt file in the source's root directory.
 *
 ********************************************************
 * Table reader:
 * reads potential file and do table splining
 *
*/

#include <fstream>
#include <sstream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cfloat>
#include <cmath>
#include "omdtool.h"
#include "param.h"
#include "treader.h"

using namespace omd;
using std::ifstream;
using std::ofstream;
using std::string;
using std::istringstream;

// Calculate spline coeficients
void TableReader::calculate_coef() {
	double m;
	int    n=n_data-2;
	
	double *b=new double[n];
	double *d=new double[n];
	double *S=new double[n_data];
	
	for (int i=1;i<n+1;i++) d[i-1]=coef[i+1].d-2.0*coef[i].d + coef[i-1].d; 
	for (int i=1;i<n-1;i++) b[i]=4.0;
	b[0]=b[n-1]=5.0;	
	
	// Calculate matrix with Thomas algorithm
	for (int k=1;k<n;k++) {m=1.0/b[k-1]; b[k]-=m; d[k]-=m*d[k-1];}
	S[n+1]=S[n]=d[n-1]/b[n-1];
	for (int k=n-2; k>=0; k--) S[k+1]=(d[k]-1.0*S[k+2])/b[k];
	S[0]=S[1];
	
	for (int i=0; i<n_data-1;i++) {
		coef[i].a=(S[i+1]-S[i])/dr/dr/dr;
		coef[i].b= 3.0*S[i]/dr/dr;
		coef[i].c=(coef[i+1].d-coef[i].d-2.0*S[i]-S[i+1])/dr;
	}
	
	delete[] b;
	delete[] d;
	delete[] S;
}

TableReader::TableReader(){
	roffset = 0.0;
	coef=NULL;
	pseu_code="#$";
	name="";
	ready=allow_outrange_low=allow_outrange_hi=false;
	outrange_hi=outrange_low=0.0;
}
                         
TableReader::TableReader(string table_filename, string table_name, 
                         double vscale, double vconst) {
	roffset = 0.0;
	coef=NULL;
	pseu_code="#$";
	ready=allow_outrange_low=allow_outrange_hi=false;
	outrange_hi=outrange_low=0.0;
	open(table_filename,table_name,vscale,vconst);
}

TableReader::TableReader(string compact_string) {
	roffset = 0.0;
	coef=NULL;
	pseu_code="#$";
	name="";
	ready=allow_outrange_low=allow_outrange_hi=false;
	outrange_hi=outrange_low=0.0;
	open(compact_string);
}

TableReader::~TableReader() {
  if(coef) delete[] coef;
}

void TableReader::allocate() {
	coef = new CoefStruct[n_data];
	for (int i=0; i<n_data; i++)
    coef[i].a=coef[i].b=coef[i].c=coef[i].d=0.0;
}

// open file and read table data
// for opening omd table file vscale&vconst shall be 1
// this is for compatibility reason

bool TableReader::open_omd(string table_filename, string table_name,
                           double vscale, double vconst) {
	
	//---------initialization-------------			
	string sformat("PLAIN");
	
	filename.assign(table_filename);
	if(table_name!="") name.assign(table_name);
	else name.assign(filename);
		
	if(!param.read_pseudo(table_filename, pseu_code+table_name)) return false;
  
	n_data=param.int_value("NumberOfData");
	allocate();
	dr=vscale*param.double_value("Spacing");   // scaling applied
	param.peek("Offset", roffset, 0.0);
	param.peek("Format", sformat);	
	
	if(param.exist("OutOfRangeLower")){
		outrange_low=param.double_value("OutOfRangeLower");
		allow_outrange_low=true;
	}

	if(param.exist("OutOfRangeUpper")){
		outrange_hi=param.double_value("OutOfRangeUpper");
		allow_outrange_hi=true;
	}
	
	if(sformat=="PLAIN") format=plain;
	else if(sformat=="RMULT") format=rmult;
	else throw string("unsuported format ")+sformat;

	ifstream fl(filename.c_str());
	if(!fl.good()) {
	  throw (string("open_omd(): error reading data file: ")+filename).c_str();
    }
	
	int lnum=0;
	char line[1024];
	bool partfound=false;
	string stag="#$"+table_name;
	
	while(fl.good()) {
		fl.getline(line, 1023); lnum++;
		istringstream sstr(line); 
		string s, ss; 
		sstr>>s>>ss;
		if(s==stag&&ss=="--") {partfound=true;break;}
	}
	if(!partfound) 
    throw string("can not find table ")+name+"@"+filename;
	
	for(int i=0;i<n_data;i++){
		fl>>coef[i].d;
		coef[i].d*=vconst;   // scaling with constant
	}

	fl.close();	
	calculate_coef();
	return true;
}

// raw table containes x y tabulated value pair
// with equispaced x.

bool TableReader::open_raw(string table_filename,double vscale, double vconst) {
  char line[2048];
  int lnum=0;
  double aold,xmin,delta=0.0;
  vector<double> vco;
  
  ifstream fl(filename.c_str());
  
  if(!fl.good()) {
    throw (string("open_raw(): error reading data file: ")+filename).c_str();
  }
  
  while(fl.good()) {
    double a,b;
    fl.getline(line,2047);
    istringstream ss(line);
    
    if(ss>>a>>b) {
      if(vco.empty()) xmin=a;
      else delta+=vscale*(a-aold); // abcissa&ordinatescaling applied
      vco.push_back(vconst*b);
      aold=a;
    }
    lnum++;
  }
  fl.close();
  if(vco.empty()) return false;
  
  // transfer
  name=table_filename.substr(0,table_filename.find('.'));
  n_data=vco.size();
  dr=delta/double(n_data-1);
  roffset=xmin;
  format=plain;
  
  allocate();
  for(int i=0;i<n_data;i++) coef[i].d=vco[i];
  calculate_coef();
  
  return true;
}

// del candidate
void TableReader::open(string table_filename, string tablename,
                       double vscale, double vconst) {
  
	if(!open_omd(table_filename,tablename,vscale,vconst)) {
		if(tablename!="") throw "the parameter insists omd format!";
		if(!open_raw(table_filename,vscale,vconst))
			throw "failed loading table";
	}
	ready=true;
}

// The format of compact string:
//     table@filename:const
//     table@filename:scale:const
//     
//     ':' may be replaced by space
// table, scale, const are optional

void TableReader::open(string compact_string) {
  istringstream ss(replace_char(compact_string,':',' '));
  string fname,tname;
  double xscl, yscl;
  ss>>fname;
  if(!(ss>>xscl)) xscl=1.0;
  if(!(ss>>yscl)) {yscl=xscl;xscl=1.0;}
  istringstream sn(replace_char(fname,'@',' '));
  sn>>tname;
  name=tname;
  if(!(sn>>fname)) {
	  fname=tname;
	  tname="";
  }
  open(fname,tname,xscl,yscl);
}

void TableReader::read(double r, double& val, double& dval) {
	
	double rval=r-roffset;
	int    x = (int)(rval/dr);
	double dx = rval-(double)x*dr;
	
	if (x>n_data-1) {
		if(allow_outrange_hi) {val=dval=outrange_hi;return;}
		else {
			die("high limit exeeded in "+name+"@"+filename+" (read) "
				"index="+as_string(x)+" r="+as_string(r)+" delta="+as_string(dx)+
				" low_limit="+as_string(min_range())+" hi_limit="+as_string(max_range()));
		}
	}
	
	if (x<-1) {
		if(allow_outrange_low) {
			val=dval=outrange_low;
			return;
		} else {
			die("low limit exeeded in "+name+"@"+filename+" (read) "
				"index="+as_string(x)+" r="+as_string(r)+" delta="+as_string(dx)+
				" low_limit="+as_string(min_range())+" hi_limit="+as_string(max_range()));
				
		}
	}
	
	if (x==n_data-1) {dx+=dr; x--;}
	if (x==-1) {x=0;}
	val=(coef[x].d + dx*(coef[x].c + dx*(coef[x].b + dx*coef[x].a)));
	dval=(coef[x].c + dx*(2.0*coef[x].b + 3.0*dx*coef[x].a));
  	
	if(format==rmult){
		val/=r;
		dval=(dval-val)/r;
	}
}

double TableReader::read(double r) {
	
	double rval=r-roffset;
	int x = (int)(rval/dr);
	double dx = rval-(double)x*dr;

	if (x>n_data-1) {
		if(allow_outrange_hi) {
			return outrange_hi;
		} else {
			die("high limit exeeded in "+name+"@"+filename+" (read) "
				"index="+as_string(x)+" r="+as_string(r)+" delta="+as_string(dx)+
 				" low_limit="+as_string(min_range())+" hi_limit="+as_string(max_range()));

		}
	}
	
	if (x<(-1)) {
		if(allow_outrange_low) {
			return outrange_low;
		} else {
			die("low limit exeeded in "+name+"@"+filename+" (read) "
				"index="+as_string(x)+" r="+as_string(r)+" delta="+as_string(dx)+
				" low_limit="+as_string(min_range())+" hi_limit="+as_string(max_range()));
		}
	}

	if (x==n_data-1) {dx+=dr; x--;}
	if (x==-1) {x=0;}
	double val=(coef[x].d + dx*(coef[x].c + dx*(coef[x].b + dx*coef[x].a)));
	if(format==rmult) val/=r;
	
	return val;
}

double TableReader::dread(double r) {
	double val, dval;
	read(r, val, dval);
	return dval;
}

double TableReader::dread2(double r) {
  
	double rval=r-roffset;
	int x = (int)(rval/dr);
	double dx = rval-(double)x*dr;
	

	if (x>n_data-1) {
		if(allow_outrange_hi) {
			return outrange_hi;
		} else {
			die("high limit exeeded in "+name+"@"+filename+" (dread2)"
				"index="+as_string(x)+" r="+as_string(r)+" delta="+as_string(dx)+
				" low_limit="+as_string(min_range())+" hi_limit="+as_string(max_range()));

		}
	}
	
	if (x<-1) {
		if(allow_outrange_low) {
			return outrange_low;
		} else {
			die("low limit exeeded in "+name+"@"+filename+" "
				"index="+as_string(x)+" r="+as_string(r)+" delta="+as_string(dx)+
				" low_limit="+as_string(min_range())+" hi_limit="+as_string(max_range()));
		}
	}

	if(x==n_data-1){dx+=dr; x--;}
	if(x==-1)x=0;

	double ddval=2.0*coef[x].b + 6.0*dx*coef[x].a;
	
	if(format==rmult) {
		double dval=dread(r);
		ddval=(ddval-2.0*dval)/r;
	}
  
	return ddval;
}

void TableReader::dump(std::ostream& ofl, int resolution) { 
	double rmax=roffset+(double)(n_data-1)*dr;
	double dx=rmax/(double)resolution;
  std::cerr<<" max="<<rmax<<" dr="<<dr<<" dx="<<dx<<std::endl;
	ofl << "# Table of value first and second derivative\n";
	for (double x=roffset; x<=rmax; x+=dx) {
		ofl << x << ' ' << read(x) << ' ' << dread(x) << ' ' << dread2(x) << '\n';
	}
}

void TableReader::dump(string filename, int resolution) { 
	ofstream ofl(filename.c_str());
  if(!ofl.good())
    throw string("can not open file for writing ")+filename;
  dump(ofl,resolution);
  ofl.close();
}

void TableReader::dump_var() {
	std::cout << "offset = " << roffset << ", " << "dr = " << dr << "\n";
}

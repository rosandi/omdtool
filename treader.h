/* omdlib
 ********************************************************
 *
 * ObjectMD
 * Molecular Dynamics Class Library
 * Version 3.0 (2009,2011) 
 *
 * Released under:
 *       GNU GENERAL PUBLIC LICENSE Version 3
 *       Copyright (C) 2005  Yudi Rosandi
 *
 * Read license.txt file in the source's root directory.
 *
 ********************************************************
 * ObjectMD header file
 *
 *      Table reader
 *      The spline data is stored in memory
 *
*/

#ifndef _TABLE_READER_
#define _TABLE_READER_
#include "param.h"

namespace omd {

  class TableReader {
    struct CoefStruct {double a, b, c, d;} *coef;
    enum tabel_format {plain, rmult} format;
    int	n_data;
    bool allow_outrange_low;
    bool allow_outrange_hi;
    double outrange_low;
    double outrange_hi;
    double  dr, roffset;
    void calculate_coef();
    void allocate();
    bool open_omd(string, string, double, double);
    bool open_raw(string, double, double);  
    bool ready;

  public:
    ParamHandler param;
    string  filename;
    string  name; //this table name may be replaced after loading file
    string  pseu_code;
    TableReader();
    TableReader(string table_filename, string table_name, double vscale=1.0, double vconst=1.0);
    TableReader(string compact_string);
    ~TableReader();
    void   open(string);
    void   open(string table_filename, string table_name, double vscale=1.0, double vconst=1.0);
    double max_range() {return (roffset+(double)(n_data-1)*dr);}
    double min_range() {return roffset;}
    bool   is_ready() {return ready;}
    void   rename(string newname){name=newname;}
    double read(double r);
    void   read(double r,double& val,double& dval);
    double dread(double r);
    double dread2(double r);
    void   dump(string filename, int resolution=1000);
    void   dump(std::ostream& ofl, int resolution=1000);
    void   dump_var();
    void   outrange(double lo, double hi) {
      outrange_low=lo;
      outrange_hi=hi;    
      allow_outrange_low=allow_outrange_hi=true; 
    }
                          
  };
  
}

#endif

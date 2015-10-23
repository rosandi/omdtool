/* drift class */

#ifndef __OMD_BUFFER_H__
#define __OMD_BUFFER_H__

#include <string>
#include <vector>
#include <fstream>

using std::string;
using std::vector;
using std::ofstream;

namespace omd {

  class Buffer {
  protected:
    bool ref; // only refference
    double *prox,*proy,*proz;
    // extra field for dump. Now osed only for vtr
    vector<string> extname;
    vector<string> vecname;
    vector<string> normname;
    
    // the duty of caller to ensure sufficient buffer length 
    vector<Buffer*> extfield;
    vector<Buffer*> vecX;
    vector<Buffer*> vecY;
    vector<Buffer*> vecZ;
    vector<Buffer*> normX;
    vector<Buffer*> normY;
    vector<Buffer*> normZ;
    
  public:
    string  name;
    bool empty; // this var name is ambigue
    double* ptr;
    int nx,ny,nz;
    
    Buffer();
    Buffer(string vname, int dx, int dy, int dz);
    void alloc(int dx, int dy, int dz);
    
    void refto(double*,int,int,int);
    void refto(Buffer* bref);
    virtual ~Buffer();
    
    int size();
    double& value(int ix, int iy, int iz);
    void dump_nz(std::ostream& ofl, int slab);
    void dump_ny(std::ostream& ofl, int slab);
    void dump_nx(std::ostream& ofl, int slab);
    void dump(std::ostream& ofl, int slab, const char norm);
    void dump(std::ostream& ofl, const char norm);
    void dump(std::ostream& ofl);
    void dumpxyz(std::ostream& ofl, string delim=" "); // xyz sparse
    void dumpcsv(std::ostream& ofl); // comma separated
    void dumpvtr(std::ostream& ofl); // vtk rectilinear
    void dumpxyz(string fname,string delim=" ");
    void dumpcsv(string fname);
    void dumpvtr(string fname);
    
    void include(string, Buffer&);
    void include_vector(string, Buffer&, Buffer&, Buffer&);
    void include_normal(string, Buffer&, Buffer&, Buffer&);
    void clear_extrafields();
    double* prof_z();
    double* prof_y();
    double* prof_x();
    
    /**
     Calculate lateral average (a profile) in axis direction.
     caller is responsible to free (using delete[]) the returned pointer
     */    
    
    double* prof(const char axis) {
      if(axis=='x') return prof_x();
      if(axis=='y') return prof_y();
      if(axis=='z') return prof_z();
      return NULL;
    }

    Buffer& add(Buffer& A);
    Buffer& add(double v);
    Buffer& sub(Buffer& A);  
    Buffer& sub(double v);
    Buffer& mul(Buffer& A);
    Buffer& mul(double v);
    Buffer& div(Buffer& A);
    Buffer& div(double v);
    Buffer& copy(Buffer& A);
    Buffer& fill(double v);
    Buffer& fill(string sym, double delta=1.0, int slab=0, int nslab=-1); 
  };
  
}

#endif

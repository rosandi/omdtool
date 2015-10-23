/* omdtool buffer class */
/* (c) Yudi Rosandi yudi@rosandi.com */

#include <iomanip>
#include <fstream>
#include <cstring>
#include <cmath>
#include "buffer.h"

using namespace omd;
using std::ofstream;

Buffer::Buffer() {
  ptr=prox=proy=proz=NULL;
  ref=true;
  empty=true;
}
  
Buffer::Buffer(string vname, int dx, int dy, int dz) {
  name=vname;
  nx=dx;ny=dy;nz=dz;
  ptr=prox=proy=proz=NULL;
  ref=false;
  alloc(dx,dy,dz);  
}

void Buffer::alloc(int dx, int dy, int dz) {
  if(dx*dy*dz) {
    ptr=new double[dx*dy*dz];
    for(int i=0;i<dx*dy*dz;i++) ptr[i]=0.0;
    empty=false;
  } else {
    ptr=NULL;
    empty=true;
  }    
}

void Buffer::refto(double* p, int _nx, int _ny, int _nz) {
  ref=true;
  nx=_nx;
  ny=_ny;
  nz=_nz;
  ptr=p;
  name.assign("noname");
  ref=true;
  empty=false;
}

void Buffer::refto(Buffer* bref) {
  nx=bref->nx;
  ny=bref->ny;
  nz=bref->nz;
  ptr=bref->ptr;
  name.assign("*");
  name.append(bref->name);
  ref=true;
  empty=false;
}

Buffer::~Buffer(){
  if(!ref) {
    if(ptr) delete[] ptr;
    if(prox) delete[] prox;
    if(proy) delete[] proy;
    if(proz) delete[] proz;
  }
}

int Buffer::size(){return nx*ny*nz;}

double& Buffer::value(int ix, int iy, int iz) {
  if((unsigned)ix>=nx || (unsigned)iy>=ny || (unsigned)iz>=nz) throw "out of range";
  return ptr[(iz*nx*ny)+ix+iy*nx];
}

void Buffer::dump_nz(std::ostream& ofl, int slab) { // default only one slab
  if(slab>=nz) throw "slab out of range";
  for(int j=0;j<ny;j++) {
    for(int i=0;i<nx;i++) 
      ofl<<value(i,j,slab)<<" ";
    ofl<<std::endl;
  }
}

void Buffer::dump_ny(std::ostream& ofl, int slab) { // default only one slab
  if(slab>=ny) throw "slab out of range";
  for(int j=0;j<nz;j++) {
    for(int i=0;i<nx;i++) 
      ofl<<value(i,slab,j)<<" ";
    ofl<<std::endl;
  }
}

void Buffer::dump_nx(std::ostream& ofl, int slab) { // default only one slab
  if(slab>=nx) throw "slab out of range";
  for(int j=0;j<nz;j++) {
    for(int i=0;i<ny;i++) 
      ofl<<value(slab,i,j)<<" ";
    ofl<<std::endl;
  }
}

void Buffer::dump(std::ostream& ofl, int slab, const char norm) {
  if(norm=='z') dump_nz(ofl,slab);
  if(norm=='y') dump_ny(ofl,slab);
  if(norm=='x') dump_nx(ofl,slab);    
}

void Buffer::dump(std::ostream& ofl, const char axis) {
  double* buf=prof(axis);
  int nn=(axis=='z')?nz:(axis=='y')?ny:nx;
  for(int i=0;i<nn;i++) {
    ofl<<i<<" "<<buf[i]<<std::endl;
  }
}

void Buffer::dump(std::ostream& ofl) {
  for(int k=0;k<nz;k++) {
    for(int j=0;j<ny;j++) {
      for(int i=0;i<nx;i++) 
        ofl<<value(i,j,k)<<" ";
      ofl<<std::endl;
    }
  }
}

void Buffer::dumpxyz(std::ostream& ofl, string delim){
  for(int k=0;k<nz;k++) {
    for(int j=0;j<ny;j++) {
      for(int i=0;i<nx;i++) {
        ofl<<i<<delim<<j<<delim<<k<<delim<<value(i,j,k);
        
        for(int v=0;v<(int)extfield.size();v++) 
          ofl<<delim<<extfield[v]->value(i,j,k);
        
        for(int v=0;v<(int)vecX.size();v++)
          ofl<<delim<<vecX[v]->value(i,j,k)
             <<delim<<vecY[v]->value(i,j,k)
             <<delim<<vecZ[v]->value(i,j,k);

        for(int v=0;v<(int)normX.size();v++) {
          double vx=normX[v]->value(i,j,k);
          double vy=normY[v]->value(i,j,k);
          double vz=normZ[v]->value(i,j,k);
          double mag=sqrt(vx*vx+vy*vy+vz*vz);
          ofl<<delim<<vx/mag<<delim<<vy/mag<<delim<<vz/mag;
        }
        
        ofl<<std::endl;
      }
    }
  }
}

void Buffer::dumpcsv(std::ostream& ofl){Buffer::dumpxyz(", ");}

void Buffer::dumpvtr(std::ostream& ofl){
  // print header
  ofl<<"<VTKFile type='RectilinearGrid' version='0.1' format='ascii'>\n"
     <<"<RectilinearGrid WholeExtent='1 "<<nx<<" 1 "<<ny<<" 1 "<<nz<<"'>\n"
     <<"<Piece Extent='1 "<<nx<<" 1 "<<ny<<" 1 "<<nz<<"'>\n"
     <<"<Coordinates>\n";
  ofl<<"<DataArray type='Float32' Name='X_COORDINATES' NumberOfComponents='1' format='ascii'>\n";
  ofl<<std::fixed<<std::setprecision(10);
  for(int i=0;i<nx;i++)ofl<<i<<" ";ofl<<"\n";
  ofl<<"</DataArray>\n";
  ofl<<"<DataArray type='Float32' Name='Y_COORDINATES' NumberOfComponents='1' format='ascii'>\n";
  ofl<<std::fixed<<std::setprecision(10);
  for(int i=0;i<ny;i++)ofl<<i<<" ";ofl<<"\n";
  ofl<<"</DataArray>\n";
  ofl<<"<DataArray type='Float32' Name='Z_COORDINATES' NumberOfComponents='1' format='ascii'>\n";
  ofl<<std::fixed<<std::setprecision(10);
  for(int i=0;i<nz;i++)ofl<<i<<" ";ofl<<"\n";
  ofl<<"</DataArray>\n</Coordinates>\n";
  ofl<<"<PointData>\n";
  
  ofl<<"<DataArray type='Float32' Name='"<<name<<"' NumberOfComponents='1' format='ascii'>\n";
  
  for(int k=0;k<nz;k++)
    for(int j=0;j<ny;j++)
      for(int i=0;i<nx;i++)
        ofl<<value(i,j,k)<<"\n";

  ofl<<"</DataArray>\n";
  
  // extra field if exist
  for(int v=0;v<(int)extfield.size();v++) {
    ofl<<"<DataArray type='Float32' Name='"<<extname[v]
       <<"' NumberOfComponents='1' format='ascii'>\n";
    
    for(int k=0;k<nz;k++)
      for(int j=0;j<ny;j++)
        for(int i=0;i<nx;i++)
          ofl<<extfield[v]->value(i,j,k)<<"\n";
    
    ofl<<"</DataArray>\n";
  }

  // extra vector field if exist
  for(int v=0;v<(int)vecX.size();v++) {
    ofl<<"<DataArray type='Float32' Name='"<<vecname[v]
       <<"' NumberOfComponents='3' format='ascii'>\n";
    
    for(int k=0;k<nz;k++)
      for(int j=0;j<ny;j++)
        for(int i=0;i<nx;i++)
          ofl<<vecX[v]->value(i,j,k)<<" "
             <<vecY[v]->value(i,j,k)<<" "
             <<vecZ[v]->value(i,j,k)<<"\n";
    
    ofl<<"</DataArray>\n";
  }

  // extra vector field if exist
  for(int v=0;v<(int)normX.size();v++) {
    ofl<<"<DataArray type='Float32' Name='"<<normname[v]
       <<"' NumberOfComponents='3' format='ascii'>\n";
    
    for(int k=0;k<nz;k++)
      for(int j=0;j<ny;j++)
        for(int i=0;i<nx;i++) {
          double vx=normX[v]->value(i,j,k);
          double vy=normY[v]->value(i,j,k);
          double vz=normZ[v]->value(i,j,k);
          double mag=sqrt(vx*vx+vy*vy+vz*vz);
          ofl<<vx/mag<<" "<<vy/mag<<" "<<vz/mag<<"\n";
        }
    ofl<<"</DataArray>\n";
  }
  
  ofl<<"</PointData>\n</Piece>\n</RectilinearGrid>\n</VTKFile>\n";
}

void Buffer::dumpxyz(string fname, string delim) {
  ofstream of(fname.c_str());
  dumpxyz(of,delim);
  of.close();
}

void Buffer::dumpcsv(string fname) {
  ofstream of(fname.c_str());
  dumpcsv(of);
  of.close();
}
void Buffer::dumpvtr(string fname) {
  ofstream of(fname.c_str());
  dumpvtr(of);
  of.close();
}

void Buffer::include(string name, Buffer& B) {
  extname.push_back(name);
  extfield.push_back(&B);
}

void Buffer::include_vector(string name, Buffer& VX, Buffer& VY, Buffer& VZ) {
  vecname.push_back(name);
  vecX.push_back(&VX);
  vecY.push_back(&VY);
  vecZ.push_back(&VZ);  
}

void Buffer::include_normal(string name, Buffer& VX, Buffer& VY, Buffer& VZ) {
  normname.push_back(name);
  normX.push_back(&VX);
  normY.push_back(&VY);
  normZ.push_back(&VZ);  
}

void Buffer::clear_extrafields() {
  extname.clear();
  vecname.clear();
  extfield.clear();
  vecX.clear();
  vecY.clear();
  vecZ.clear();
}

Buffer& Buffer::add(Buffer& A) {
  if(nx!=A.nx || ny!=A.ny) throw "incompatible size";
  for(int i=0;i<size();i++) {
    ptr[i]+=A.ptr[i];
  }
  return *this;
}

Buffer& Buffer::add(double v) {
  for(int i=0;i<size();i++) ptr[i]+=v;
  return *this;
}

Buffer& Buffer::sub(Buffer& A) {
  if(nx!=A.nx || ny!=A.ny) throw "incompatible size";
  for(int i=0;i<size();i++) {
    ptr[i]-=A.ptr[i];
  }
  return *this;
}

Buffer& Buffer::sub(double v) {
  for(int i=0;i<size();i++) ptr[i]-=v;
  return *this;
}

Buffer& Buffer::mul(Buffer& A) {
  if(nx!=A.nx || ny!=A.ny) throw "incompatible size";
  for(int i=0;i<size();i++) {
    ptr[i]*=A.ptr[i];
  }
  return *this;
}  

Buffer& Buffer::mul(double v) {
  for(int i=0;i<size();i++) ptr[i]*=v;
  return *this;
}  

Buffer& Buffer::div(Buffer& A) {
  if(nx!=A.nx || ny!=A.ny) throw "incompatible size";
  for(int i=0;i<size();i++) {
    if(A.ptr[i]==0.0) throw "division by zero";
    ptr[i]/=A.ptr[i];
  }
  return *this;
}

Buffer& Buffer::div(double v) {
  if(v==0.0) throw "division by zero";
  for(int i=0;i<size();i++) ptr[i]/=v;
  return *this;
}  

Buffer& Buffer::copy(Buffer& A) {
  if(nx!=A.nx || ny!=A.ny) throw "incompatible size";
  for(int i=0;i<size();i++) {
    ptr[i]=A.ptr[i];
  }
  return *this;
}

Buffer& Buffer::fill(double v) {
  int sz=size();
  for(int i=0;i<sz;i++) ptr[i]=v;
  return *this;
}

Buffer& Buffer::fill(string sym, double delta, int slab, int nslab) {
  if(nslab<0) nslab=nz-slab;
  else nslab=slab+nslab;
  if(nslab>nz) throw "slab out of range";
  
  if(sym=="x") {
    for(int k=slab;k<nslab;k++) {
      for(int j=0;j<ny;j++)
        for(int i=0;i<nx;i++)
          value(i,j,k)=(double)i*delta;
    }
  }
  
  else if(sym=="y") {
    for(int k=0;k<nz;k++) {
      for(int j=0;j<ny;j++)
        for(int i=0;i<nx;i++)
          value(i,j,k)=(double)j*delta;
    }
  }
  
  else if(sym=="z") {
    for(int k=0;k<nz;k++) {
      for(int j=0;j<ny;j++)
        for(int i=0;i<nx;i++)
          value(i,j,k)=(double)k*delta;
    }
  }
  
  else throw "unknown symbol";
  
  return *this;
}

double*  Buffer::prof_z() {
  if(!proz) proz=new double[nz];
  memset(proz,0,nz*sizeof(double));
  double nc=nx*ny;
  
  for(int k=0;k<nz;k++) {
    for(int j=0;j<ny;j++) 
      for(int i=0;i<nx;i++)
        proz[k]+=value(i,j,k);
    proz[k]/=nc;
  }
  
  return proz;
}

double*  Buffer::prof_y() {
  if(!proy) proy=new double[ny];
  memset(proy,0,ny*sizeof(double));
  double nc=nx*nz;
  
  for(int j=0;j<ny;j++) {
    for(int k=0;k<nz;k++) 
      for(int i=0;i<nx;i++)
        proy[j]+=value(i,j,k);
    proy[j]/=nc;
  }
  
  return proy;
}

double*  Buffer::prof_x() {
  if(!prox) prox=new double[nx];
  memset(prox,0,nx*sizeof(double));
  double nc=nz*ny;
  
  for(int i=0;i<nx;i++) {
    for(int j=0;j<ny;j++) 
      for(int k=0;k<nz;k++)
        prox[i]+=value(i,j,k);
    prox[i]/=nc;
  }
  
  return prox;
}

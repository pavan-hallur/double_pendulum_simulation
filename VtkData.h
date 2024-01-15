#ifndef vtkHeader_H
#define vtkHeader_H

#include <vector>
#include <string>

struct Vec3
{
	double x;
	double y;
	double z;

	Vec3() : x(0), y(0), z(0) {}
	Vec3(const double x, const double y, const double z) : x(x), y(y), z(z) {}
};

class VtkData
{
private:

    std::vector< Vec3 > points;
    
    std::vector< std::vector< int > > cells;
    
    std::vector< double > scalarPointData;
    std::vector< Vec3 >   vectorPointData;
    
    bool printWarnings = true;
    
public:
    
    void setPoints ( const std::vector< Vec3 > points );
    
    void setCells ( const std::vector< std::vector< int > > cells );
    
    void setScalarPointData ( const std::vector< double > scalarPointData );
    void setVectorPointData ( const std::vector< Vec3 >   vectorPointData );
    
    void writeVtkFile ( const std::string filename, const bool verbose = true ) const;
    
    void setPrintWarnings ( const bool printWarnings );
};

#endif
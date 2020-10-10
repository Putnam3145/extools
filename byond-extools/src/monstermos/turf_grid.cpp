#include "turf_grid.h"
#include "../dmdism/opcodes.h"
#include "var_ids.h"
#include <algorithm>
#include <cstring>
#include <memory>
#include <execution>
#include <complex>
using namespace monstermos::constants;

extern Value SSair;

template<class T>
void add_sources(Eigen::SparseMatrix<T>& x, Eigen::SparseMatrix<T>& x0)
{
	x += x0;
}

/// technically unstable; just don't give it a rate of more than 1/4
template<class T>
void diffuse(Eigen::SparseMatrix<T>& g, Eigen::SparseMatrix<T>& g0,ProcessSpecifier& processing_turfs,float rate = 0.1)
{
	for(int y=0;y<ATMOS_WIDTH;++y)
	{
		for(int x=0;x<ATMOS_WIDTH;++x)
		{
			if(processing_turfs(y,x) > ProcessingLevel::PROCESSING_LEVEL_SPACE)
			{
				T coeff;
				float shared = 0.0;
				if(x>0 && processing_turfs(y,x-1))
				{
					coeff += g0.coeff(y,x-1);
					shared++;
				}
				if(x<254 && processing_turfs(y,x+1))
				{
					coeff += g0.coeff(y,x+1);
					shared++;
				}
				if(y>0 && processing_turfs(y-1,x))
				{
					coeff += g0.coeff(y-1,x);
					shared++;
				}
				if(y<254 && processing_turfs(y+1,x))
				{
					coeff += g0.coeff(y+1,x);
					shared++;
				}
				g.coeffRef(y,x) = g0.coeff(y,x)+(coeff - shared*g0.coeff(y,x))*rate;
			}
		}
	}
}

template<class T>
void advect(Eigen::SparseMatrix<T>& g,Eigen::SparseMatrix<T>& g0,VelocityMatrix& v,ProcessSpecifier& processing_turfs)
{
	for(int y=0;y<ATMOS_WIDTH;y++)
	{
		for(int x=0;x<ATMOS_WIDTH;x++)
		{
			if(processing_turfs.coeff(y,x))
			{
				std::complex<float> d = v.coeff(y,x);
				if(d != 0.0f)
				{
					float nx = x-d.real();
					float ny = y-d.imag();
					nx = std::clamp(nx,0.0f,static_cast<float>(ATMOS_WIDTH));
					ny = std::clamp(nx,0.0f,static_cast<float>(ATMOS_WIDTH));
					float s1 = nx - (int)nx;
					float s0 = 1 - s1;
					float t1 = ny - (int)ny;
					float t0 = 1 - t1;
					g.coeffRef(y,x) = s0*(t0*g0.coeff(ny,nx))+t1*g0.coeff(ny+1,nx)+
									s1*(t0*g0.coeff(ny,nx+1))+t1*g0.coeff(ny+1,nx+1);
				}
			}
		}
	}
}

void project(VelocityMatrix& v,ProcessSpecifier& processing_turfs)
{
	Eigen::ArrayXf p(ATMOS_WIDTH,ATMOS_WIDTH);
	Eigen::ArrayXf div(ATMOS_WIDTH,ATMOS_WIDTH);
	p.setZero();
	float h = 1.0f/ATMOS_WIDTH;
	for(int k=0;k<v.outerSize();++k)
	{
		for(VelocityMatrix::InnerIterator it(v,k);it;++it)
		{
			auto x = it.col();
			auto y = it.row();
			div(x,y) = -0.5*(v.coeff(x+1,y).real()-v.coeff(x-1,y).real()) + 
							(v.coeff(x,y+1).imag()-v.coeff(x,y-1).imag());
		}
	}
	for(int k=0;k<20;k++)
	{
		for(int i=0;i<ATMOS_WIDTH;i++)
		{
			for(int j=0;j<ATMOS_WIDTH;j++)
			{
				p(i,j) = (div(i,j)+p(i-1,j)+p(i+1,j)+p(i,j+1)+p(i,j-1))/4;
			}
		}
		p = (processing_turfs < 2).select(0.0f,p);
	}
	for(int y=0;y<ATMOS_WIDTH;y++)
	{
		for(int x=0;x<ATMOS_WIDTH;x++)
		{
			if(processing_turfs(y,x) > ProcessingLevel::PROCESSING_LEVEL_SPACE)
			{
				using namespace std::complex_literals;
				auto& ref = v.coeffRef(y,x);
				ref -= (0.5f*((p(y,x-1)-p(y,x+1))+(p(y+1,x)-p(y-1,x)*1if)));
			}
		}
	}
}

void TurfGrid::process() {
	gasesByTurf.clear();
	for(int i=0;i<maxz;i++)
	{
		auto& proc = processing_turfs.at(i);
		VelocityMatrix& v = gasVelocity[i];
		for(auto& p : gasDensities)
		{
				GasMatrix& g = gasDensities[p.first][i];
				GasMatrix& g0 = sources[p.first][i];
				add_sources(g,g0);
				diffuse(g0,g,proc);
				advect(g,g0,v,proc);
		}
		{
			GasMatrix& e = gasEnergy[i];
			GasMatrix& e0 = energySources[i];
			add_sources(e,e0);
			diffuse(e0,e,proc);
			advect(e,e0,v,proc);
		}
		{
			VelocityMatrix& v0 = velocitySources[i];
			add_sources(v,v0);
			diffuse(v0,v,proc);
			project(v0,proc);
			advect(v,v0,v,proc);
			project(v,proc);
		}
	}
}

GasMixture TurfGrid::gas_from_turf(int x,int y,int z,bool cache)
{
	return gas_from_turf(x+y*maxx+z*maxx*maxy,cache);
}

GasMixture TurfGrid::gas_from_turf(int id,bool cache)
{
	if(gasesByTurf.contains(id))
	{
		return gasesByTurf.at(id);
	}
	GasMixture newMix(2500);
	int x = id%maxx;
	int y = (id/maxx)%maxy;
	int z = (id/(maxy*maxx));
	for(auto& g : gasDensities)
	{
		auto actualGas = g.second[z].coeff(y,x);
		if(actualGas > 0.0f)
			newMix.set_moles(g.first,actualGas);
	}
	newMix.set_temperature(gasEnergy[z].coeff(y,x)/newMix.heat_capacity());
	if(cache)
	{
		gasesByTurf.insert({id,newMix});
	}
	return newMix;
}

void TurfGrid::refresh() {
	maxx = Value::World().get("maxx").valuef;
	maxy = Value::World().get("maxy").valuef;
	maxz = Value::World().get("maxz").valuef;
	maxid = maxx * maxy * maxz;
	gasEnergy.resize(maxz);
	gasVelocity.resize(maxz);
	velocitySources.resize(maxz);
	energySources.resize(maxz);
	processing_turfs.resize(maxz);
	for(int i=0;i<maxz;i++)
	{
		for(auto& p : gasDensities)
		{
			p.second[i].conservativeResize(ATMOS_WIDTH,ATMOS_WIDTH);
			sources[p.first][i].conservativeResize(ATMOS_WIDTH,ATMOS_WIDTH);
		}
		gasEnergy[i].conservativeResize(ATMOS_WIDTH,ATMOS_WIDTH);
		gasVelocity[i].conservativeResize(ATMOS_WIDTH,ATMOS_WIDTH);
		velocitySources[i].conservativeResize(ATMOS_WIDTH,ATMOS_WIDTH);
		energySources[i].conservativeResize(ATMOS_WIDTH,ATMOS_WIDTH);
		processing_turfs[i].conservativeResize(ATMOS_WIDTH,ATMOS_WIDTH);
	}
}
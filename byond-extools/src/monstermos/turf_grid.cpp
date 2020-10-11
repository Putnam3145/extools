#include "turf_grid.h"
#include "../dmdism/opcodes.h"
#include "var_ids.h"
#include <algorithm>
#include <cstring>
#include <memory>
#include <execution>
#include <complex>
#include <compare>
#define _USE_MATH_DEFINES
#include <math.h>
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
	for(int k=0;k<processing_turfs.outerSize();++k)
	{
		for(ProcessSpecifier::InnerIterator it(processing_turfs,k);it;++it)
		{
			if(it.value() > ProcessingLevel::TURF_PROCESSING_LEVEL_SPACE)
			{
				int x = it.col();
				int y = it.row();
				T coeff;
				float shared = 0.0;
				if(x>0 && processing_turfs(y,x-1) > ProcessingLevel::TURF_PROCESSING_LEVEL_SPACE)
				{
					coeff += g0.coeff(y,x-1);
					shared++;
				}
				if(x<254 && processing_turfs(y,x+1) > ProcessingLevel::TURF_PROCESSING_LEVEL_SPACE)
				{
					coeff += g0.coeff(y,x+1);
					shared++;
				}
				if(y>0 && processing_turfs(y-1,x) > ProcessingLevel::TURF_PROCESSING_LEVEL_SPACE)
				{
					coeff += g0.coeff(y-1,x);
					shared++;
				}
				if(y<254 && processing_turfs(y+1,x) > ProcessingLevel::TURF_PROCESSING_LEVEL_SPACE)
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
	for(int k=0;k<processing_turfs.outerSize();++k)
	{
		for(ProcessSpecifier::InnerIterator it(processing_turfs,k);it;++it)
		{
			int x = it.col();
			int y = it.row();
			std::complex<float> d = v.coeff(y,x);
			if(d != 0.0f)
			{
				float nx = x-d.real();
				float ny = y-d.imag();
				nx = std::clamp(nx,0.0f,static_cast<float>(g.cols()));
				ny = std::clamp(nx,0.0f,static_cast<float>(g.rows()));
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

void project(VelocityMatrix& v,ProcessSpecifier& processing_turfs)
{
	Eigen::ArrayXf p;
	Eigen::ArrayXf div;
	p.resizeLike(v);
	div.resizeLike(v);
	float h = 1.0f/v.rows();
	for(int k=0;k<v.outerSize();++k)
	{
		for(VelocityMatrix::InnerIterator it(v,k);it;++it)
		{
			auto x = it.col();
			auto y = it.row();
			div(y,x) = -0.5*(v.coeff(y,x+1).real()-v.coeff(y,x-1).real()) + 
							(v.coeff(y+1,x).imag()-v.coeff(y-1,x).imag());
		}
	}
	for(int k=0;k<20;k++)
	{
		for(int j=0;j<processing_turfs.outerSize();++k)
		{
			for(ProcessSpecifier::InnerIterator it(processing_turfs,j);it;++it)
			{
				if(it.value() > ProcessingLevel::TURF_PROCESSING_LEVEL_SPACE)
				{
					int x = it.col();
					int y = it.row();
					p(y,x) = (div(y,x)+p(y-1,x)+p(y+1,x)+p(y,x+1)+p(y,x-1))/4;
				}
			}
		}
		p = (processing_turfs > 0).select(0.0f,p);
	}
	for(int k=0;k<processing_turfs.outerSize();++k)
	{
		for(ProcessSpecifier::InnerIterator it(processing_turfs,k);it;++it)
		{
			if(it.value() > ProcessingLevel::TURF_PROCESSING_LEVEL_SPACE)
			{
				int x = it.col();
				int y = it.row();
				using namespace std::complex_literals;
				auto& ref = v.coeffRef(y,x);
				ref -= (0.5f*((p(y,x-1)-p(y,x+1))+(p(y+1,x)-p(y-1,x)*1if)));
			}
		}
	}
}

void stretch(VelocityMatrix& v,ProcessSpecifier& processing_turfs)
{
	bool done = false;
	while(!done)
	{
		done = true;
		for(int k=0;k<v.outerSize();++k)
		{
			for(VelocityMatrix::InnerIterator it(v,k);it;++it)
			{
				if(std::abs(it.value()) > MAX_ATMOS_VELOCITY)
				{
					auto original_value = it.value();
					it.valueRef() /= std::abs(it.value());
					const auto velocityLeftover = original_value - it.value();
					if(std::abs(velocityLeftover) > MAX_ATMOS_VELOCITY) done = false;
					const auto normed = it.value(); // already normalized it earlier
					if(velocityLeftover.real() < 0)
					{
						v.coeffRef(it.col()-1,it.row()) += (velocityLeftover.real())*(normed.real()*normed.real());
					}
					else if(velocityLeftover.real() > 0)
					{
						v.coeffRef(it.col()+1,it.row()) += (velocityLeftover.real())*(normed.real()*normed.real());
					}
					if(velocityLeftover.imag() < 0)
					{
						v.coeffRef(it.col(),it.row()-1) += (velocityLeftover.imag())*(normed.imag()*normed.imag());
					}
					else if(velocityLeftover.imag() > 0)
					{
						v.coeffRef(it.col(),it.row()+1) += (velocityLeftover.imag())*(normed.imag()*normed.imag());
					}
				}
				else if(processing_turfs.coeff(it.col(),it.row()) == 0)
				{
					it.valueRef() = std::complex(0.0f);
				}
			}
		}
	}
}

void ZAtmos::resize(int y,int x) {
	gasEnergy.conservativeResize(y,x);
	gasVelocity.conservativeResize(y,x);
	velocitySources.conservativeResize(y,x);
	energySources.conservativeResize(y,x);
	processing_turfs.conservativeResize(y,x);
}

void TurfGrid::process() {
	gasesByTurf.clear();
	std::for_each(
		std::execution::par_unseq,
		zLevels.begin(),zLevels.end(),
		[](ZAtmos& level)
	{
		ProcessSpecifier& proc = level.processing_turfs;
		VelocityMatrix& v = level.gasVelocity;
		{
			VelocityMatrix& v0 = level.velocitySources;
			add_sources(v,v0);
			stretch(v,proc);
			diffuse(v0,v,proc,0.005);
			project(v0,proc);
			advect(v,v0,v,proc);
			project(v,proc);
			v0.setZero();
		}
		std::for_each(std::execution::par_unseq,level.sources.begin(),level.sources.end(),[&](auto& p)
		{
			GasMatrix& g = level.gasDensities[p.first];
			GasMatrix& g0 = p.second;
			add_sources(g,g0);
			diffuse(g0,g,proc,0.05);
			advect(g,g0,v,proc);
			g0.setZero();
		});
		{
			Eigen::SparseMatrix<float>& e = level.gasEnergy;
			Eigen::SparseMatrix<float>& e0 = level.energySources;
			add_sources(e,e0);
			diffuse(e0,e,proc,0.1);
			advect(e,e0,v,proc);
			e0.setZero();
		}
	});
}

VelocityMatrix getGradient(const GasMatrix& p,ProcessSpecifier& boundary)
{
	std::vector<Eigen::Triplet<std::complex<float>>> res;
	res.reserve(boundary.nonZeros());
	for(int k=0;k<p.outerSize();++k)
	{
		for(GasMatrix::InnerIterator it(p,k);it;++it)
		{
			auto x = it.col();
			auto y = it.row();
			if(boundary(y,x) > TURF_PROCESSING_LEVEL_SPACE)
			{
				auto base = it.value();
				auto nx = (base-p.coeff(y,x+1))-(base-(p.coeff(y,x-1)));
				auto ny = (base-p.coeff(y+1,x))-(base-(p.coeff(y-1,x)));
				res.emplace_back(y,x,std::complex<float>(nx,ny));
			}
		}
	}
	VelocityMatrix trueRes;
	trueRes.setFromTriplets(res.begin(),res.end());
	return trueRes;
}

void TurfGrid::consider_pressure_differences()
{
	std::for_each(
		std::execution::par_unseq,
		zLevels.begin(),zLevels.end(),
		[](ZAtmos& level) {
		ProcessSpecifier& proc = level.processing_turfs;
		GasMatrix heatCapacities;
		GasMatrix totals;
		for(auto& p : level.gasDensities)
		{
			GasMatrix& g = p.second;
			heatCapacities += g * gas_specific_heat[p.first];
			totals += g;
		}
		GasMatrix temperatures = level.gasEnergy.cwiseProduct(heatCapacities);
		GasMatrix pressures = ((temperatures.cwiseProduct(totals)) * R_IDEAL_GAS_EQUATION) / (CELL_VOLUME/(ATMOS_GRANULARITY*ATMOS_GRANULARITY));
		VelocityMatrix pressureGradients = getGradient(pressures,proc);
		pressureGradients /= (ONE_ATMOSPHERE);
		level.velocitySources += pressureGradients;
	});
}

GasMixture TurfGrid::gas_from_turf(int x,int y,int z,bool cache)
{
	return gas_from_turf(x+y*maxx+z*maxx*maxy,cache);
}

template<class T,int size=ATMOS_GRANULARITY>
T getBlockSum(Eigen::SparseMatrix<T>& mat,int y,int x)
{
	return mat.block(y,x,size,size).sum();
}

GasMixture TurfGrid::gas_from_turf(int id,bool cache)
{
	if(gasesByTurf.contains(id))
	{
		return gasesByTurf.at(id);
	}
	GasMixture newMix(2500);
	int x = (id%maxx)*ATMOS_GRANULARITY;
	int y = ((id/maxx)%maxy)*ATMOS_GRANULARITY;
	int z = (id/(maxy*maxx));
	auto& zLevel = zLevels[z];
	for(auto& g : zLevel.sources) // sources is guaranteed to have all of 'em, even roundstart
	{
		float gasAmount = getBlockSum(g.second,y,x) + getBlockSum(zLevel.gasDensities[g.first],y,x);
		if(gasAmount > 0.0f)
			newMix.set_moles(g.first,gasAmount);
	}
	newMix.set_temperature((getBlockSum(zLevel.gasEnergy,y,x)+getBlockSum(zLevel.energySources,y,x))/newMix.heat_capacity());
	if(cache)
	{
		gasesByTurf.insert({id,newMix});
	}
	return newMix;
}

void TurfGrid::add_to_source(int id,GasMixture& gas,float force,float dir)
{
	int x = (id%maxx)*ATMOS_GRANULARITY;
	int y = ((id/maxx)%maxy)*ATMOS_GRANULARITY;
	int z = (id/(maxy*maxx));
	auto& zLevel = zLevels[z];
	auto ener = gas.thermal_energy();
	float theta = 0.0;
	for(auto i = 0;i<ATMOS_GRANULARITY;i++)
	{
		for(auto j=0;j<ATMOS_GRANULARITY;j++)
		{
			for(auto& g : gas.gases_in_mix())
			{
				zLevel.sources.try_emplace(g,maxy*ATMOS_GRANULARITY,maxx*ATMOS_GRANULARITY);
				zLevel.sources[g].coeffRef(y+j,x+i) += gas.get_moles(g);
			}
			zLevel.energySources.coeffRef(y+j,x+i) += ener;
			zLevel.velocitySources.coeffRef(y+j,x+i) += std::polar(force,dir);
		}
	}
}

void TurfGrid::set_turf_proc_level(int id, int level)
{
	int x = (id%maxx)*ATMOS_GRANULARITY;
	int y = ((id/maxx)%maxy)*ATMOS_GRANULARITY;
	int z = (id/(maxy*maxx));
	auto& zLevel = zLevels[z];
	for(auto i = 0;i<ATMOS_GRANULARITY;i++)
	{
		for(auto j=0;j<ATMOS_GRANULARITY;j++)
		{
			zLevel.processing_turfs.coeffRef(y+j,x+i) = level;
		}
	}
}

std::complex<float> TurfGrid::get_velocity(int id)
{
	int x = (id%maxx)*ATMOS_GRANULARITY;
	int y = ((id/maxx)%maxy)*ATMOS_GRANULARITY;
	int z = (id/(maxy*maxx));
	auto& zLevel = zLevels[z];
	std::complex<float> ret;
	for(auto i = 0;i<ATMOS_GRANULARITY;i++)
	{
		for(auto j=0;j<ATMOS_GRANULARITY;j++)
		{
			ret += zLevel.gasVelocity.coeff(y+j,x+i);
		}
	}
	return ret;
}

void TurfGrid::refresh() {
	maxx = Value::World().get("maxx").valuef;
	maxy = Value::World().get("maxy").valuef;
	maxz = Value::World().get("maxz").valuef;
	maxid = maxx * maxy * maxz;
	zLevels.resize(maxz);
	std::for_each(zLevels.begin(),zLevels.end(),[&](auto& zLevel){
		zLevel.resize(maxy*ATMOS_GRANULARITY,maxx*ATMOS_GRANULARITY);
	});
}
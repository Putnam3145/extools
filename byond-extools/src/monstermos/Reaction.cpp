#include "Reaction.h"

#include <algorithm>

using namespace monstermos::constants;

Reaction::Reaction(Value v)
{
    List min_reqs = v.get("min_requirements");
    if(min_reqs.at("TEMP").type == DataType::NUMBER) min_temp_req = min_reqs.at("TEMP");
    if(min_reqs.at("ENER").type == DataType::NUMBER) min_ener_req = min_reqs.at("ENER");
    if(min_reqs.at("MAX_TEMP").type == DataType::NUMBER) max_temp_req = min_reqs.at("MAX_TEMP");
    for(unsigned int i=0;i < total_num_gases;i++)
    {
        auto gasReq = min_reqs.at(gas_id_to_type[i]);
        if(gasReq.type == DataType::NUMBER)
        {
            min_gas_reqs[i] = gasReq;
        }
    }
    major_gas = v.get("major_gas");
    priority = v.get("priority");
    auto proc = Core::try_get_proc(Core::stringify(v.get("type")) + "/react");
    if(!proc)
    {
        Core::alert_dd("Could not find proc for reaction! " + Core::stringify(v.get("type")) + "/react");
    }
    else
    {
        Core::alert_dd("Found proc for reaction: " + Core::stringify(v.get("type")) + "/react");
    }
    proc_id = proc->id;
}

bool Reaction::check_conditions(const GasMixture& air) const
{
    return (air.get_moles(major_gas) > 0.0 &&
            air.get_temperature() >= min_temp_req && air.get_temperature() <= max_temp_req && 
            min_ener_req > 0.0 && air.thermal_energy() >= min_ener_req &&
            std::all_of(min_gas_reqs.cbegin(),min_gas_reqs.cend(),[&](auto& info) {
                return air.get_moles(info.first) >= info.second;
            }));
}

int Reaction::react(GasMixture& air,Value src,Value holder) const
{
    return (int)(float)(Core::get_proc(proc_id).call({src,holder}));
}
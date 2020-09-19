#define NOMINMAX
#include <iostream>
#include <functional>
#include <algorithm>

#include "disassembler.h"
#include "opcodes.h"
#include "helpers.h"
#include "context.h"


Disassembler::Disassembler(std::vector<std::uint32_t> bc, const std::vector<Core::Proc>& ps)
{
	context_ = std::make_unique<Context>(bc, ps);
}

Disassembler::Disassembler(std::uint32_t* bc, unsigned int bc_len, const std::vector<Core::Proc>& ps)
{
	std::vector<std::uint32_t> v(bc, bc + bc_len);
	context_ = std::make_unique<Context>(v, ps);
}

Disassembly Disassembler::disassemble()
{
	std::vector<Instruction> instrs;
	instrs.reserve(512);
	while (context_->more()) {
		instrs.push_back(disassemble_next());
	}

	return Disassembly(std::move(instrs));
}

Instruction Disassembler::disassemble_next()
{
	auto offset = context_->current_offset();
	auto root = context_->eat(nullptr);
	Instruction instr { (Bytecode) root };

	instr.set_offset(offset);
	instr.add_byte(root);

	DisassembleCallback callback = get_disassemble_callback(root);
	if (callback)
	{
		callback(&instr, context(), this);
	}

	return instr;
}

bool Disassembler::disassemble_var(Instruction& instr)
{
	switch ((AccessModifier) context_->peek())
	{
	case AccessModifier::SUBVAR:
	{
		std::uint32_t val = context_->eat(&instr);
		instr.opcode().add_info(" SUBVAR");
		if (!disassemble_var(instr))
		{
			return true;
		}
		instr.add_comment(".");
		if (!disassemble_var(instr))
		{
			return true;
		}

		break;
	}

	case AccessModifier::LOCAL:
	case AccessModifier::GLOBAL:
	//case CACHE:
	case AccessModifier::ARG:
	{
		std::uint32_t type = context_->eat(&instr);
		std::uint32_t localno = context_->eat(&instr);

		std::string modifier_name = "UNKNOWN_MODIFIER";
		if (modifier_names.find(static_cast<AccessModifier>(type)) != modifier_names.end())
		{
			modifier_name = modifier_names.at(static_cast<AccessModifier>(type));
		}

		instr.opcode().add_info(" " + modifier_name + std::to_string(localno));
		instr.add_comment(modifier_name + std::to_string(localno));
		break;
	}
	case AccessModifier::INITIAL:
	{
		context_->eat(&instr);
		instr.add_comment("INITIAL(");
		std::uint32_t val = context_->eat(&instr);
		instr.add_comment(byond_tostring(val));
		instr.add_comment(")");
		break;
	}
	case AccessModifier::CACHE:
	{
		context_->eat(&instr);
		instr.add_comment("CACHE");
		break;
	}
	case AccessModifier::WORLD:
	case AccessModifier::NULL_:
	case AccessModifier::DOT:
	case AccessModifier::SRC:
	{
		std::uint32_t type = context_->eat(&instr);

		std::string modifier_name = "UNKNOWN_MODIFIER";
		if (modifier_names.find(static_cast<AccessModifier>(type)) != modifier_names.end())
		{
			modifier_name = modifier_names.at(static_cast<AccessModifier>(type));
		}

		instr.opcode().add_info(" " + modifier_name);
		instr.add_comment(modifier_name);
		break;
	}
	case AccessModifier::ARGS:
		context_->eat(&instr);
		instr.add_comment("ARGS");
		break;
	case AccessModifier::PROC_NO_RET:
	case AccessModifier::PROC:
	case AccessModifier::SRC_PROC:
	case AccessModifier::SRC_PROC_SPEC:
		return false;
	default:
	{
		std::uint32_t val = context_->eat(&instr);
		instr.add_comment(byond_tostring(val));
		break;
	}
	}
	return true;
}

bool Disassembler::disassemble_proc(Instruction& instr)
{
	switch ((AccessModifier) context_->peek())
	{
		case AccessModifier::PROC_NO_RET:
		case AccessModifier::PROC:
		{
			context_->eat(&instr);
			instr.add_comment(Core::get_proc(context_->eat(&instr)).simple_name);
			add_call_args(instr, context_->eat(&instr));
			return true;
		}
		case AccessModifier::SRC_PROC:
		case AccessModifier::SRC_PROC_SPEC:
		{
			context_->eat(&instr);
			std::string name = byond_tostring(context_->eat(&instr));
			std::replace(name.begin(), name.end(), ' ', '_');
			instr.add_comment(name);
			add_call_args(instr, context_->eat(&instr));
			return true;
		}
		default:
			instr.add_comment("<!>");
			return false;
	}
}

void Disassembler::add_call_args(Instruction& instr, unsigned int num_args)
{
	instr.add_comment(" [num_args: " + std::to_string(num_args) + "]");
}

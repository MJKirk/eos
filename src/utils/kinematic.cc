/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Danny van Dyk
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <src/utils/kinematic.hh>
#include <src/utils/private_implementation_pattern-impl.hh>
#include <src/utils/stringify.hh>

#include <map>
#include <vector>

namespace eos
{
    template <>
    struct Implementation<Kinematics>
    {
        std::vector<double> variables_data;

        std::map<std::string, unsigned> variables_map;
    };

    Kinematics::Kinematics() :
        PrivateImplementationPattern<Kinematics>(new Implementation<Kinematics>)
    {
    }

    Kinematics::~Kinematics()
    {
    }

    Kinematics
    Kinematics::clone() const
    {
        Kinematics result;

        result._imp.reset(new Implementation<Kinematics>(*_imp));

        return result;
    }

    KinematicVariable
    Kinematics::operator[] (const std::string & name) const
    {
        auto i(_imp->variables_map.find(name));

        if (_imp->variables_map.end() == i)
            throw UnknownKinematicVariableError(name);

        return KinematicVariable(_imp, i->second);
    }

    std::string
    Kinematics::as_string () const
    {
        std::string result;

        auto i(_imp->variables_map.cbegin()), i_end(_imp->variables_map.cend());
        if (i != i_end)
        {
            result = i->first + '=' + stringify(_imp->variables_data[i->second]);
            ++i;
        }

        for ( ; i != i_end ; ++i)
        {
            result += ", " + i->first + '=' + stringify(_imp->variables_data[i->second]);
        }

        return result;
    }

    KinematicVariable
    Kinematics::declare(const std::string & name, const double & value)
    {
        auto i(_imp->variables_map.find(name));

        if (_imp->variables_map.end() == i)
        {
            int index = _imp->variables_data.size();
            _imp->variables_map[name] = index;
            _imp->variables_data.push_back(value);

            return KinematicVariable(_imp, index);
        }

        _imp->variables_data[i->second] = value;

        return KinematicVariable(_imp, i->second);
    }

    void
    Kinematics::set(const std::string & name, const double & value)
    {
        auto i(_imp->variables_map.find(name));

        if (_imp->variables_map.end() == i)
            throw UnknownKinematicVariableError(name);

        _imp->variables_data[i->second] = value;
    }

    KinematicVariable::KinematicVariable(const std::shared_ptr<Implementation<Kinematics>> & imp, unsigned index) :
        _imp(imp),
        _index(index)
    {
    }

    KinematicVariable::~KinematicVariable()
    {
    }

    KinematicVariable::operator double () const
    {
        return _imp->variables_data[_index];
    }

    double
    KinematicVariable::operator() () const
    {
        return _imp->variables_data[_index];
    }

    const KinematicVariable &
    KinematicVariable::operator= (const double & value)
    {
        _imp->variables_data[_index] = value;

        return *this;
    }

    UnknownKinematicVariableError::UnknownKinematicVariableError(const std::string & variable) throw () :
        Exception("Unknown kinematic variable: '" + variable + "'")
    {
    }
}

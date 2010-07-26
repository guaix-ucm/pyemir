/*
 * Copyright 2008-2010 Sergio Pascual
 *
 * This file is part of PyEmir
 *
 * PyEmir is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PyEmir is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PyEmir.  If not, see <http://www.gnu.org/licenses/>.
 *
 */


#include "reject_factory.h"
#include "reject_methods.h"

namespace Numina {

auto_ptr<RejectMethod>
RejectMethodFactory::create(const std::string& name,
		PyObject* args,
		auto_ptr<CombineMethod> combine_method) {
	if (name == "none")
		return auto_ptr<RejectMethod>(new NoneReject(args, combine_method));
	return auto_ptr<RejectMethod>();
}


} // namespace Numina


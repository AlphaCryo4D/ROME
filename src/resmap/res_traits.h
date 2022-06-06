/***************************************************************************
 *
 * Intel® Parallel Computing Center for Structural Biology
 * Principal Investigator : Youdong (Jack) Mao (Youdong_Mao@dfci.harvard.edu)
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
 *
 * Authors: "Jian Wang(wj_hust08@hust.edu.cn)"
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/

#pragma once

#define JN_SELF(T) (*(T*)(nullptr))

namespace cppsci {

    template<bool A, bool B> struct static_or                 { enum { value = true  }; };
    template<>               struct static_or<false, false>   { enum { value = false }; };
#define JN_STATIC_OR(A, B)   static_or<A, B>::value

    template<bool A, bool B> struct static_and                { enum { value = false }; };
    template<>               struct static_and<true, true>    { enum { value = true  }; }; 
#define JN_STATIC_AND(A, B)  static_and<A, B>::value

    template<bool A>         struct static_not                { enum { value = false }; };
    template<>               struct static_not<false>         { enum { value = true  }; }; 
#define JN_STATIC_NOT(A)     static_not<A>::value

} // namespace cppsci


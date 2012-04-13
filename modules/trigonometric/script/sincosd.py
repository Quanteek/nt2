[ ## this file was manually modified by jt
    {
     'functor' : {
         'description' :['simultaneous sine and cosine from angle in degree',
                         '\\\\arg tie(si,co) = sincosd(a) returns the sine and cosine of a.',
                         '\\\\arg si = sincosd(a,co) returns the sine of a and put the cosine of a in co.',
                         '\\\\arg sincosd(a,si,co) returns 0 and and put the sine of a in si and the cosine cosine of a in co.',
                         ],
         'param_0' : ['angle in degree'],   
          'special_synopsis' : [
                "template<class A0> inline",
                "typename boost::dispatch::meta::call<tag::sincosd(A0 const&,A0 const&)",
                "                                    >::type",
                "sincosd(A0 const& a);",
                " ",
                "template<class A0> inline",
                "A0 sincosd(A0 const& a,A0,A0& si);",
                " ",
                "template<class A0> inline",
                "int sincosd(A0 const& a,A0,A0& si,A0& co);",
                 ],
         'arity' : '1',
         'call_types' : [],
         'ret_arity' : '2',
         'rturn' : {
             'default' : 'boost::fusion::tuple<ftype,ftype>',
            },
         'simd_types' : ['real_'],
         'special' : ['trigonometric'],
         'type_defs' : ['typedef typename boost::dispatch::meta::as_floating<T>::type ftype;'],
         'types' : ['real_', 'unsigned_int_', 'signed_int_'],
        },
     'info' : 'manually modified',
     'unit' : {
         'global_header' : {
             'first_stamp' : 'created  by jt the 11/02/2011',
             'included' : 
                ['#include <nt2/toolbox/trigonometric/constants.hpp>',
                 '#include <nt2/include/functions/sind.hpp>',
                 '#include <nt2/include/functions/cosd.hpp>'],
             'notes' : [],
             'stamp' : 'modified by jt the 11/02/2011',
            },
         'ranges' : {
             'default' : [['T(-2000)', 'T(2000)']],
            },
         'specific_values' : {
             'default' : {
                 'nt2::Zero<T>()' : {'result' : ['nt2::Zero<r_t0>()', 'nt2::One<r_t0>()'],'ulp_thresh' : '0.75',},
                 'nt2::_45<T>()' : {'result' : ['nt2::Sqrt_2o_2<r_t0>()', 'nt2::Sqrt_2o_2<r_t0>()'],'ulp_thresh' : '0.75',},
                 'nt2::_90<T>()' : {'result' : ['nt2::One<r_t0>()', 'nt2::Zero<r_t0>()'],'ulp_thresh' : '0.75',},
                },
             'real_' : {
                 '-nt2::_180<T>()' : {'result' : ['nt2::Zero<r_t0>()', 'nt2::Mone<r_t1>()'],'ulp_thresh' : '0.75',},
                 '-nt2::_45<T>()' : {'result' : ['-nt2::Sqrt_2o_2<r_t0>()', 'nt2::Sqrt_2o_2<r_t0>()'],'ulp_thresh' : '0.75',},
                 '-nt2::_90<T>()' : {'result' : ['nt2::Mone<r_t0>()', 'nt2::Zero<r_t1>()'],'ulp_thresh' : '0.75',},
                 'nt2::Inf<T>()' : {'result' : ['nt2::Nan<r_t0>()', 'nt2::Nan<r_t0>()'],'ulp_thresh' : '0.75',},
                 'nt2::Minf<T>()' : {'result' : ['nt2::Nan<r_t0>()', 'nt2::Nan<r_t0>()'],'ulp_thresh' : '0.75',},
                 'nt2::Nan<T>()' : {'result' : ['nt2::Nan<r_t0>()', 'nt2::Nan<r_t0>()'],'ulp_thresh' : '0.75',},
                 'nt2::Zero<T>()' : {'result' : ['nt2::Zero<r_t0>()', 'nt2::One<r_t0>()'],'ulp_thresh' : '0.75',},
                 'nt2::_180<T>()' : {'result' : ['nt2::Zero<r_t0>()', 'nt2::Mone<r_t1>()'],'ulp_thresh' : '0.75',},
                 'nt2::_45<T>()' : {'result' : ['nt2::Sqrt_2o_2<r_t0>()', 'nt2::Sqrt_2o_2<r_t0>()'],'ulp_thresh' : '0.75',},
                 'nt2::_90<T>()' : {'result' : ['nt2::One<r_t0>()', 'nt2::Zero<r_t0>()'],'ulp_thresh' : '0.75',},
                },
             'signed_int_' : {
                 '-nt2::_45<T>()' : {'result' : ['-nt2::Sqrt_2o_2<r_t0>()', 'nt2::Sqrt_2o_2<r_t0>()'],'ulp_thresh' : '0.75',},
                 '-nt2::_90<T>()' : {'result' : ['nt2::Mone<r_t0>()', 'nt2::Zero<r_t1>()'],'ulp_thresh' : '0.75',},
                 'nt2::Zero<T>()' : {'result' : ['nt2::Zero<r_t0>()', 'nt2::One<r_t0>()'],'ulp_thresh' : '0.75',},
                 'nt2::_45<T>()' : {'result' : ['nt2::Sqrt_2o_2<r_t0>()', 'nt2::Sqrt_2o_2<r_t0>()'],'ulp_thresh' : '0.75',},
                 'nt2::_90<T>()' : {'result' : ['nt2::One<r_t0>()', 'nt2::Zero<r_t0>()'],'ulp_thresh' : '0.75',},
                },
            },
         'verif_test' : {
             'property_call' : {
                 'real_' : ['nt2::sincosd(a0)'],
                },
             'property_value' : {
                 'real_' : [['nt2::sind(a0)', 'nt2::cosd(a0)']],
                },
             'ulp_thresh' : {
                 'real_' : ['0.5'],
                },
            },
        },
    },
]

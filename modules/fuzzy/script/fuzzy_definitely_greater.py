[ ## this file was manually modified by jt
    {
     'functor' : {
            'description' : ['\\\\f$a_0 \\\\gt a_1+a_2*\max(|a_0|,|a1|)\\\\f$'],
         'arity' : '3',
         'call_types' : [],
         'ret_arity' : '0',
         'rturn' : {
             'default' : 'typename nt2::meta::as_logical<T>::type',
            },
         'simd_types' : ['real_', 'signed_int_', 'unsigned_int_'],
         'special' : ['fuzzy'],
         'type_defs' : [],
         'types' : ['real_', 'signed_int_', 'unsigned_int_'],
         'simd_types' : ['real_'],         
        },
     'info' : 'manually modified',
     'unit' : {
         'global_header' : {
             'first_stamp' : 'created  by jt the 04/03/2011',
             'included' : ['#include <nt2/sdk/simd/logical.hpp>'],
             'no_ulp' : 'True',
             'notes' : [],
             'stamp' : 'modified by jt the 04/03/2011',
            },
         'ranges' : {
             'default' : [['T(-10)', 'T(10)'], ['T(-10)', 'T(10)'], ['T(-10)', 'T(10)']],
             'unsigned_int_' : [['T(0)', 'T(10)'], ['T(0)', 'T(10)'], ['iT(0)', 'iT(10)']],
            },
         'specific_values' : {
             'default' : {
                 'T(0),T(0),T(1)' : {'result' : 'nt2::False<r_t>()','ulp_thresh' : '0',},
                 'T(0),T(1),T(1)' : {'result' : 'nt2::False<r_t>()','ulp_thresh' : '0',},
                },
            },
         'verif_test' : {
             'property_call' : {
                 'default' : ['nt2::fuzzy_definitely_greater(a0,a1,a2)'],
                },
             'property_value' : {
                 'default' : ['nt2::fuzzy_definitely_greater(a0,a1,a2)'],
                },
             'simd' : {
                },
             'ulp_thresh' : {
                 'default' : ['0'],
                },
            },
        },
    },
]

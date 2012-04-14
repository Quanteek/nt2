/// these functions are lacking and prevent other to perform or
/// even to be compiled properly

///                            | Factorisations            | general             | solvers                   | gallery,  clement
/// diag               used in | plu, pqr, qr, schur, svd  | trace               | solve_svd_ip              |
/// expand             used in | plu, pqr, qr, schur, svd  |                     |                           |
/// prod               used in | plu                       |                     |                           |
/// sum                used in | plu, pqr, qr              | norm, trace         |                           |
/// _(i, j) indexation used in | plu, pqr, qr, svd         |                     | solve_qr_ip, solve_svd_ip | cauchy, chebspec
/// _(i, _) indexation used in |                           |                     |                           | chebvand
/// allto              used in |                           |                     |                           | cycol
/// allfrom            used in |                           |                     |                           | dorr
/// trans              used in | plu, svd                  | rot90               |                           | magic
/// * (matrix mul)     used in | svd                       | normest, subspace   |                           |
/// maximum (max)      used in |                           | norm                |                           |
/// minimum (min)      used in |                           | vecnorm             |                           |
/// norm_eucl          used in |                           | normest             |                           |
/// rand               used in |                           | normest             |                           |
/// flipud             used in |                           | rot90               |                           | circul
/// normp              used in |                           | vecnorm             |                           |
/// bsxfun             used in |                           | vecnorm             |                           | cauchy, fielder
/// cath               used in |                           |                     |                           | cycol
/// catv               used in |                           |                     |                           | magic
/// repmat             used in |                           |                     |                           | vandermonde

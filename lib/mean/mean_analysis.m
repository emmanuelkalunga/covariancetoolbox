function [M,sum_dis] = mean_analysis(COVMAT, alpha, method)

switch method
    case 'riemann'
        [M,sum_dis] = riemann_mean_analysis(COVMAT);
    case 'sdivergence'
        [M,sum_dis] = sdiv_mean_analysis(COVMAT);
    case 'ld'
        [M,sum_dis] = logdet_mean_alpha_analysis(COVMAT,alpha);
    case 'bhat'
        [M,sum_dis] = logdet_mean_alpha_analysis(COVMAT,alpha);
    case 'opttransp'
        [M,sum_dis] = opttransp_mean_analysis(COVMAT);
end
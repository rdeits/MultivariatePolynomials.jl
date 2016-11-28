using MultivariatePolynomials
using BenchmarkTools

let
    @polyvar x y
    p = x + y + 2x^2 + y^3

    vars = [x, y]
    vals = [1.0, 2.0]

    println(@benchmark(($p)($vals, $vars)))
    vals = MultivariatePolynomials.evalmap(MultivariatePolynomials.vars(p), vals, vars)
    @code_warntype(MultivariatePolynomials.termeval(p[1], vals))
    @code_warntype(zero(p))
    @code_warntype(p(vals, vars))
end

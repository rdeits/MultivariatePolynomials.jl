# I do not use it but I import the function to add a method
export differentiate

"""
    differentiate(p::AbstractPolynomialLike, v::AbstractVariable, deg::Int=1)

Differentiate `deg` times the polynomial `p` by the variable `v`.

To differentiate by multiple variables or with multiple degrees simultaneously, you can simply use Julia's built-in function broadcasting:

### Examples

```julia
p = 3x^2*y + x + 2y + 1
differentiate(p, x) # should return 6xy + 1

# Using the `.()` notation for broadcast:
differentiate.(p, (x, y)) # should return (6xy+1, 3x^2+1)
differentiate.(p, [x, y], 2) # should return [6y, 0]
differentiate.(p, [x, y], (1:2)') # should return: [6xy+1 6y; 3x^2+1 0]
```
"""
function differentiate end

function differentiate(p, x::Tuple, deg::Integer)
    Base.depwarn(
"""Calling `differentiate()` with a tuple of arguments is deprecated.
Instead, just use broadcasting:

To differentiate p w.r.t x and y once, returning (dp/dx, dp/dy)

    differentiate.(p, (x, y))

To differentiate p w.r.t. x and y twice, returning [dp/dx d2p/dx2; dp/dy d2p/dy2]:

    differentiate.(p, (x, y), (1:2)')
""", :MP_differentiate_tuple)
    differentiate.(p, x, (1:deg)')
end

function differentiate(p, x::AbstractVector, deg::Integer)
    Base.depwarn(
"""Calling `differentiate()` with a vector of arguments is deprecated.
Instead, just use broadcasting:

To differentiate p w.r.t x and y once, returning [dp/dx, dp/dy]

    differentiate.(p, [x, y])

To differentiate p w.r.t. x and y twice, returning [dp/dx d2p/dx2; dp/dy d2p/dy2]:

    differentiate.(p, [x, y], (1:2)')
""", :MP_differentiate_vector)
    differentiate.(p, x, (1:deg)')
end


# # Fallback for everything else
# _diff_promote_op(::Type{T}, ::Type{<:AbstractVariable}) where T = T
differentiate(α::T, v::AbstractVariable) where T = zero(T)

# _diff_promote_op(::Type{<:AbstractVariable}, ::Type{<:AbstractVariable}) = Int
differentiate(v1::AbstractVariable, v2::AbstractVariable) = v1 == v2 ? 1 : 0

# _diff_promote_op(::Type{TT}, ::Type{<:AbstractVariable}) where {T, TT<:AbstractTermLike{T}} = changecoefficienttype(TT, Base.promote_op(*, T, Int))
differentiate(t::AbstractTermLike, v::AbstractVariable) = coefficient(t) * differentiate(monomial(t), v)

# _diff_promote_op(::Type{PT}, ::Type{<:AbstractVariable}) where {T, PT<:APL{T}} = polynomialtype(PT, Base.promote_op(*, T, Int))
# The polynomial function will take care of removing the zeros
differentiate(p::APL, v::AbstractVariable) = polynomial(differentiate.(terms(p), v), SortedState())

differentiate(p::RationalPoly, v::AbstractVariable) = (differentiate(p.num, v) * p.den - p.num * differentiate(p.den, v)) / p.den^2

const ARPL = Union{APL, RationalPoly}

# _vec_diff_promote_op(::Type{PT}, ::AbstractVector{VT}) where {PT, VT} = _diff_promote_op(PT, VT)
# _vec_diff_promote_op(::Type{PT}, ::NTuple{N, VT}) where {PT, N, VT}   = _diff_promote_op(PT, VT)
# _vec_diff_promote_op(::Type{PT}, ::VT, xs...) where {PT, VT}          = _diff_promote_op(PT, VT)
# _vec_diff_promote_op(::Type{PT}, xs::Tuple) where PT = _vec_diff_promote_op(PT, xs...)

# even if I annotate with ::Array{_diff_promote_op(T, PolyVar{C}), N+1}, it cannot detect the type since it seems to be unable to determine the dimension N+1 :(
# function differentiate(ps::AbstractArray{PT, N}, xs::Union{AbstractArray, Tuple}) where {N, PT<:ARPL}
#     qs = Array{_vec_diff_promote_op(PT, xs), N+1}(uninitialized, length(xs), size(ps)...)
#     cartesian = CartesianIndices(ps)
#     for (i, x) in enumerate(xs)
#         for j in linearindices(ps)
#             J = cartesian[j]
#             qs[i, J] = differentiate(ps[J], x)
#         end
#     end
#     qs
# end

function differentiate(p::ARPL, xs)
    Base.depwarn(
"""Calling `differentiate()` with a vector of arguments is deprecated.
Instead, just use broadcasting:

To differentiate p w.r.t. x and y once, returning [dp/dx, dp/dy]

    differentiate.(p, [x, y])

To differentiate p w.r.t. x and y twice, returning [dp/dx d2p/dx2; dp/dy d2p/dy2]:

    differentiate.(p, [x, y], (1:2)')
    """, :MP_differentiate_ARPL)
   [differentiate(p, x) for x in xs]
end

function differentiate(ps::AbstractVector, xs::Union{AbstractVector, Tuple})
    Base.depwarn(
"""Calling `differentiate()` with a vector of arguments is deprecated. Instead, just use broadcasting:

To differentiate p1 and p2 w.r.t. x and y, returning [dp1/dx dp1/dy; dp2/dx dp2/dy]:

    differentiate.([p1, p2], [x, y]')
    """, :MP_differentiate_ps_xs)
    differentiate.(ps, collect(xs)')
end

# differentiate(p, [x, y]) with TypedPolynomials promote x to a Monomial
differentiate(p::ARPL, m::AbstractMonomial) = differentiate(p, variable(m))

# # In Julia v0.5, Base.promote_op returns Any for PolyVar, Monomial and MatPolynomial
# # Even on Julia v0.6 and Polynomial, Base.promote_op returns Any...
# _diff_promote_op(::Type{PT}, ::Type{VT}) where {PT, VT} = Base.promote_op(differentiate, PT, VT)
# _diff_promote_op(::Type{MT}, ::Type{<:AbstractVariable}) where {MT<:AbstractMonomialLike} = termtype(MT, Int)

function differentiate(p, x, deg::Integer)
    # # To ensure type stability, we will always return whatever
    # # the type of differentiate(p, x) would be, even when deg=0. 
    # # We compute that type here: 
    T = typeof(differentiate(p, x))
    # # Note that, as long as inference succeeds, the Julia compiler
    # # is able to optimize out the call to `differentiate(p, x)` and 
    # # just replace T with a constant.

    if deg < 0
        throw(DomainError())
    elseif deg == 0
        convert(T, p)::T
    else
        convert(T, differentiate(differentiate(p, x), x, deg-1))::T
    end
end

# function differentiate(p, x, deg::Int)
#     if deg < 0
#         throw(DomainError(deg))
#     elseif deg == 0
#         # Need the conversion with promote_op to be type stable for PolyVar, Monomial and MatPolynomial
#         @show p typeof(p)
#         @show x typeof(x)
#         @show _diff_promote_op(typeof(p), typeof(x)) 
#         return convert(_diff_promote_op(typeof(p), typeof(x)), p)
#     else
#         return differentiate(differentiate(p, x), x, deg-1)
#     end
# end

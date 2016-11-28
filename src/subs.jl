function evalmap{T}(vars, x::Vector{T}, varorder::Vector{PolyVar})
  vals = Vector{T}(length(vars))
  for (i, var) in enumerate(varorder)
    j = findfirst(vars, var)
    # If i == 0, that means that the variable is not present
    # so it is ignored
    if j > 0
      vals[j] = x[i]
    end
  end
  vals
end

function termeval(t::Term, vals::Vector)
  @assert length(vals) > 0
  val = vals[1] ^ t.x.z[1]
  for i in 2:length(vals)
      if t.x.z[i] > 0
          val *= vals[i]^t.x.z[i]
      end
  end
  val * t.α
end

function (m::PolyVar)(x::Vector, varorder)
  Term(m)(x, varorder)
end

function (m::Monomial)(x::Vector, varorder)
  Term(m)(x, varorder)
end

function (t::Term)(x::Vector, varorder)
  vals = evalmap(vars(t), x, varorder)
  termeval(t, vals)
end

function (p::VecPolynomial)(x::Vector, varorder)
  vals = evalmap(vars(p), x, varorder)
  sum(t -> termeval(t, vals), p)
end

function (p::MatPolynomial)(x::Vector, varorder)
    VecPolynomial(p)(x, varorder)
end

function (q::RationalPoly)(x::Vector, varorder)
  q.num(x, varorder) / q.den(x, varorder)
end

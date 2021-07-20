# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    signarea(A, B, C)

Compute signed area of triangle formed
by points `A`, `B` and `C`.
"""
function signarea(A::Point{2}, B::Point{2}, C::Point{2})
  ((B - A) × (C - A)) / 2
end

"""
    area(triangle)

Compute the area of `triangle`.
"""
function area(t::Triangle{3})
  vs = vertices(t)

  norm((vs[2] - vs[1]) × (vs[3] - vs[1])) / 2
end

"""
    sideof(point, segment)

Determines on which side of the oriented `segment`
the `point` lies. Possible results are `:LEFT`,
`:RIGHT` or `:ON` the segment.
"""
function sideof(p::Point{2,T}, s::Segment{2,T}) where {T}
  a, b = vertices(s)
  area = signarea(p, a, b)
  ifelse(area > atol(T), :LEFT, ifelse(area < -atol(T), :RIGHT, :ON))
end

"""
    sideof(point, chain)

Determines on which side of the closed `chain` the
`point` lies. Possible results are `:INSIDE` or
`:OUTSIDE` the chain.
"""
function sideof(p::Point{2,T}, c::Chain{2,T}) where {T}
  w = windingnumber(p, c)
  ifelse(isapprox(w, zero(T), atol=atol(T)), :OUTSIDE, :INSIDE)
end

"""
    normal(triangle)

Determine the normalised normal of `triangle`
"""
function normal(t::Triangle)
  vs = vertices(t)

  # the normal of a Triangle is the cross product of two arbitrary edges
  (vs[2] - vs[1]) × (vs[3] - vs[1]) / norm((vs[2] - vs[1]) × (vs[3] - vs[1]))
end
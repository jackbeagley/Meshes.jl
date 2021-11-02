# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    intersecttype(s1, s2)

Compute the intersection type of two line segments `s1` and `s2`.

The intersection type can be one of five types according to
Balbes, R. and Siegel, J. 1990:

1. intersect at point interior to both segments
2. intersect at end point of one and an interior point of the other
3. intersect at end point of both segments
4. overlap at more than one point
5. do not overlap nor intersect

## References

* Balbes, R. and Siegel, J. 1990. [A robust method for calculating
  the simplicity and orientation of planar polygons]
  (https://www.sciencedirect.com/science/article/abs/pii/0167839691900198)
"""
function intersecttype(s1::Segment{2,T}, s2::Segment{2,T}) where {T}
  x̄ = centroid(s1)
  ȳ = centroid(s2)
  x1, x2 = vertices(s1)
  y1, y2 = vertices(s2)
  Q = [x1, y1, x2, y2, x1]

  # winding number of x̄
  windingx = zero(T)
  determinatex = true
  for j in 1:length(Q)-1
    α = ∠(Q[j], x̄, Q[j+1])
    if α ≈ 0 || α ≈ π || α ≈ -π
      determinatex = false
    end
    windingx += α
  end

  # winding number of ȳ
  windingy = zero(T)
  determinatey = true
  for j in 1:length(Q)-1
    β = ∠(Q[j], ȳ, Q[j+1])
    if β ≈ 0 || β ≈ π || β ≈ -π
      determinatey = false
    end
    windingy += β
  end

  # winding number ≈ 0?
  wxzero = isapprox(windingx, zero(T), atol=atol(T))
  wyzero = isapprox(windingy, zero(T), atol=atol(T))

  # there are three cases to consider
  # Balbes, R. and Siegel, J. 1990.
  if determinatex && determinatey # CASE (I)
    if !wxzero && !wyzero
      # configuration (1)
      CrossingSegments(intersectpoint(s1, s2))
    else
      # configuration (5)
      NoIntersection()
    end
  elseif determinatex || determinatey # CASE (II)
    if !(determinatex ? wxzero : wyzero)
      if x1 ≈ y1 || x1 ≈ y2 || x2 ≈ y1 || x2 ≈ y2
        # configuration (3)
        CornerTouchingSegments(intersectpoint(s1, s2))
      else
        # configuration (2)
        MidTouchingSegments(intersectpoint(s1, s2))
      end
    else
      # configuration (5)
      NoIntersection()
    end
  elseif !determinatex && !determinatey # CASE (III)
    if !isapprox((x2 - x1) × (y2 - y1), zero(T), atol=atol(T)^2)
      # configuration (3)
      CornerTouchingSegments(intersectpoint(s1, s2))
    else
      # configuration (3), (4) or (5)
      intersectcolinear(s1, s2)
    end
  end
end

"""
    intersecttype(s1, s2)

Compute the intersection type of two line segments `s1` and `s2` in 3D.
See https://stackoverflow.com/questions/55220355/how-to-detect-whether-two-segmentin-3d-spaceintersect.
"""
function intersecttype(s1::Segment{3,T}, s2::Segment{3,T}) where {T}
  # get coordinates
  p₁, p₂ = coordinates.(vertices(s1))
  q₁, q₂ = coordinates.(vertices(s2))

  # create matrices and vectors to allow for evaluation of intersection
  A = [(p₂ - p₁) (q₁ - q₂)]
  b = q₁ - p₁

  λ = A \ b

  # calculate the rank of the augmented matrix
  rₐ = rank([A b])
  # calculate the rank of the rectangular matrix
  r = rank(A)

  # use matrix rank to determine basic intersection properties
  # segments are co-planar (but not co-linear)
  if (rₐ == 2) && (r == 2)
    # if either element is approximately 0 or 1, set them as so to prevent any domain errors
    λ₁ = isapprox(λ[1], zero(T), atol=atol(T)) ? zero(T) : (isapprox(λ[1], one(T), atol=atol(T)) ? one(T) : λ[1])
    λ₂ = isapprox(λ[2], zero(T), atol=atol(T)) ? zero(T) : (isapprox(λ[2], one(T), atol=atol(T)) ? one(T) : λ[2])

    # if λs are outside of the interval [0, 1], they do not intersect
    if any((λ₁, λ₂) .< zero(T)) || any((λ₁, λ₂) .> one(T))
      NoIntersection()
    # if both λs are either 0 or 1, they are CornerTouchingSegments
    elseif ((λ₁ ≈ zero(T)) || (λ₁ ≈ one(T))) && ((λ₂ ≈ zero(T)) || (λ₂ ≈ one(T)))
      CornerTouchingSegments(s1(λ₁))
    # if either λ is 0 or 1 then they are MidTouchingSegments
    elseif ((λ₁ ≈ zero(T)) || (λ₁ ≈ one(T))) ⊻ ((λ₂ ≈ zero(T)) || (λ₂ ≈ one(T)))
      MidTouchingSegments(s1(λ₁))
    # otherwise they are simple CrossingSegments
    else
      CrossingSegments(s1(λ₁))
    end
  # segments are co-linear
  elseif (rₐ == 1) && (r == 1)
    intersectcolinear(s1, s2)
  else
    NoIntersection()
  end
end

# compute the intersection of two line segments assuming that it is a point
function intersectpoint(s1::Segment{2}, s2::Segment{2})
  x1, x2 = vertices(s1)
  y1, y2 = vertices(s2)
  intersectpoint(Line(x1, x2), Line(y1, y2))
end

# intersection of two line segments assuming that they are colinear
function intersectcolinear(s1::Segment{Dim,T}, s2::Segment{Dim,T}) where {Dim,T}
  m1, M1 = coordinates.(vertices(s1))
  m2, M2 = coordinates.(vertices(s2))

  # make sure that segment vertices are "ordered"
  m1, M1 = any(m1 .> M1) ? (M1, m1) : (m1, M1)
  m2, M2 = any(m2 .> M2) ? (M2, m2) : (m2, M2)

  # relevant vertices
  u = Point(max.(m1, m2))
  v = Point(min.(M1, M2))

  if isapprox(u, v, atol=atol(T))
    CornerTouchingSegments(u)
  elseif any(coordinates(u) .< coordinates(v))
    OverlappingSegments(Segment(u, v))
  else
    NoIntersection()
  end
end

function intersecttype(s::Segment{Dim,T}, m::Mesh{Dim,T}) where {Dim, T}
  intersection_found = false
  t_intersection = NaN
  intersecting_face = 0
  intersecting_point = Point{Dim, T}
  
  # get the faces in the mesh
  m_faces = faces(m, 2)

  # iterate through each face and check for intersections (looking for the closest one)
  for (i, face) ∈ enumerate(m_faces)
      # evaluate whether or not the ray intersects the mesh
      intersection_result = intersecttype(s, face)

      if intersection_result isa IntersectingSegmentTri
          intersection_found = true

          # if this intersection is closer than a previously found one
          # flag it as so
          if !(intersection_result.t >= t_intersection)
              intersecting_face = i
              t_intersection = intersection_result.t
              intersecting_point = intersection_result.value
          end
      end
  end

  # output the closest intersection if found, and don't otherwise
  if intersection_found
      return IntersectingSegmentMesh(intersecting_point, t_intersection, intersecting_face)
  else
      return NoIntersection()
  end    
end

function intersecttype(s::Segment{3,T}, m::Plane{3,T}) where {T}
  s_v = s.vertices
  
  n = m.v × m.w
  l_dot_n = (s_v[2] - s_v[1]) ⋅ n
  pₒ_dot_n = m.pₒ.coords ⋅ n

  if isapprox(l_dot_n, zero(T))
      if isapprox(pₒ_dot_n, l_dot_n)
          return ContainedSegmentPlane()
      else
          return NoIntersection()
      end
  else
    λ = (l_dot_n - pₒ_dot_n) / l_dot_n

      if (λ < zero(T)) || (λ > one(T))
          p = s(λ)

          return IntersectingSegmentPlane(p, λ)
      else
          return NoIntersection()
      end
  end
end

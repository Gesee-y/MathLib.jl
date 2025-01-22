## Some function for the Rects type ##

"""
	mutable struct Rect2D
		x :: Int
		y :: Int

		# The Rect dimension
		w :: Int
		h :: Int
"""
mutable struct Rect2D
	x :: Int
	y :: Int

	# The Rect dimension
	w :: Int
	h :: Int
	
	# Constructors #
	
	Rect2D() = new(0,0,0,0)
	Rect2D(x::Integer,y::Integer,w::Integer,h::Integer) = new(x,y,w, h,)
	
	function Rect2D(v1::Vector2D,v2::Vector2D)
		tl = (v1.x < v2.x) ? v1 : v2
		br = (tl == v1) ? v2 : v1
		
		new(tl.x,tl.y,br.x - tl.x, br.y - tl.y)
	end
end

function point_in_rect(p::Vector2D,r::Rect2D)
	condition1 = (r.x <= p.x <= (r.x+r.w))
	condition2 = (r.y <= p.y <= (r.y+r.h))

	return (condition1 && condition2)
end

function has_intersection(r1::Rect2D,r2::Rect2D)
	dist = iVec2(r1.x-r2.x,r1.y-r2.y)
	v1 = iVec2(2r1.x-r1.w,2r1.y-r1.h)
	v2 = iVec2(2r2.x-r2.w,2r2.y-r2.h)

	return norm(dist) < (norm(v1) + norm(v2))
end

function intersection(r1::Rect2D,r2::Rect2D)
	tl = _is_the_leftmost(r1,r2) ? iVec2(r2.x,r2.y) : iVec2(r1.x,r1.y)
	br = _is_the_leftmost(r1,r2) ? iVec2(r1.x+r1.w,r1.y+r1.h) : iVec2(r2.x+r2.w,r2.y+r2.h)

	if (tl.x > br.x || tl.y > br.y) return nothing
	elseif (tl.x == br.x || tl.y == br.y)
		if tl.x == br.x & tl.y == br.y
			return tl
		end

		return Ray{Int,2}(tl,br)
	else
		return Rect2D(tl,br)
	end
end

_is_the_leftmost(r1::Rect2D,r2::Rect2D) = r1.x < r2.x
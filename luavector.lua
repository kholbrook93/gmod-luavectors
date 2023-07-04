AddCSLuaFile()

-- LuaVector implementation
-- see https://github.com/ValveSoftware/source-sdk-2013/blob/master/sp/src/mathlib/mathlib_base.cpp#L3102

local LuaVector = {}
LuaVector.__index = LuaVector
LuaVector.__type = "LuaVector"

function LuaVector:__index (k, v)
  if k == "x" then
    return rawget(self, 1)
  elseif k == "y" then
    return rawget(self, 2)
  elseif k == "z" then
    return rawget(self, 3)
  else
    local t = rawget(self, k)
    if t then return t end
    return rawget(LuaVector, k)
  end
end

-- Constants
LuaVector.FLT_EPSILON = 1.19209290e-07
LuaVector.DBL_EPSILON = 2.2204460492503131e-16
LuaVector.NORMALIZE_ANGLE_DP = 4

--- Constructor for LuaVector objects.
-- @param x (number) The x-coordinate of the vector. Default is 0.
-- @param y (number) The y-coordinate of the vector. Default is 0.
-- @param z (number) The z-coordinate of the vector. Default is 0.
-- @return (table) The new LuaVector object.
function LuaVector.new(x, y, z)
  if isvector(x) then
    return LuaVector.fromVector(x)
  end

  x = x or 0
  y = y or 0
  z = z or 0

  assert(isnumber(x), tostring(x) .. " is not a double")
  assert(isnumber(y), tostring(y) .. " is not a double")
  assert(isnumber(z), tostring(z) .. " is not a double")

  local self = setmetatable({}, LuaVector)
  self[1] = x or 0
  self[2] = y or 0
  self[3] = z or 0
  return self
end

--- Creates a new LuaVector object from a Vector object.
-- @param vec (Vector) The Vector object to create a LuaVector from.
-- @return (table) The new LuaVector object.
function LuaVector.fromVector(vec)
  return LuaVector.new(vec[1], vec[2], vec[3])
end

-- Metamethods
function LuaVector.__add(vector1, vector2)
  return LuaVector.new(
    vector1[1] + vector2[1],
    vector1[2] + vector2[2],
    vector1[3] + vector2[3]
  )
end

function LuaVector.__sub(vector1, vector2)
  return LuaVector.new(
    vector1[1] - vector2[1],
    vector1[2] - vector2[2],
    vector1[3] - vector2[3]
  )
end

function LuaVector.__mul(vector, multiplier)
  if type(multiplier) == "number" then
    return LuaVector.new(
      vector[1] * multiplier,
      vector[2] * multiplier,
      vector[3] * multiplier
    )
  elseif type(multiplier) == "table" and multiplier[1] and multiplier[2] and multiplier[3] then
    return LuaVector.new(
      vector[1] * multiplier[1],
      vector[2] * multiplier[2],
      vector[3] * multiplier[3]
    )
  else
    error("Invalid multiplier type")
  end
end

function LuaVector.__div(vector, divisor)
  return LuaVector.new(
    vector[1] / divisor,
    vector[2] / divisor,
    vector[3] / divisor
  )
end

function LuaVector.__eq(vector1, vector2)
  return vector1:IsEqualTol(vector2)
end

function LuaVector.__tostring(vector)
  return "LuaVector(" .. vector[1] .. ", " .. vector[2] .. ", " .. vector[3] .. ")"
end

--- Adds the components of another LuaVector object to this LuaVector object.
-- @param vector (table) The LuaVector object to add.
function LuaVector:Add(vector)
  self[1] = self[1] + vector[1]
  self[2] = self[2] + vector[2]
  self[3] = self[3] + vector[3]
  self.cache = nil
end

--- Subtracts the components of another LuaVector object from this LuaVector object.
-- @param vector (table) The LuaVector object to subtract.
function LuaVector:Sub(vector)
  self[1] = self[1] - vector[1]
  self[2] = self[2] - vector[2]
  self[3] = self[3] - vector[3]
  self.cache = nil
end

--- Multiplies this LuaVector object by a scalar or another LuaVector object.
-- @param multiplier (number or table) The scalar or LuaVector object to multiply by.
function LuaVector:Mul(multiplier)
  local t = type(multiplier)
  if t == "number" then
    self[1] = self[1] * multiplier
    self[2] = self[2] * multiplier
    self[3] = self[3] * multiplier
  elseif (t == "table" or t == "Vector") and multiplier[1] and multiplier[2] and multiplier[3] then
    self[1] = self[1] * multiplier[1]
    self[2] = self[2] * multiplier[2]
    self[3] = self[3] * multiplier[3]
  else
    error("Invalid multiplier type")
  end
  self.cache = nil
end

--- Divides this LuaVector object by a divisor.
-- @param divisor (number) The divisor to divide by.
function LuaVector:Div(divisor)
  self[1] = self[1] / divisor
  self[2] = self[2] / divisor
  self[3] = self[3] / divisor
  self.cache = nil
end

--- Calculates the length of this LuaVector object.
-- @return (number) The length of the LuaVector.
function LuaVector:Length()
  return math.sqrt(self[1] * self[1] + self[2] * self[2] + self[3] * self[3])
end

--- Calculates the squared length of this LuaVector object.
-- @return (number) The squared length of the LuaVector.
function LuaVector:LengthSqr()
  return self[1] * self[1] + self[2] * self[2] + self[3] * self[3]
end

--- Normalizes this LuaVector object to have a length of 1.
function LuaVector:Normalize()
  local length = self:Length()
  if length > 0 then
    self[1] = self[1] / length
    self[2] = self[2] / length
    self[3] = self[3] / length
  end
  self.cache = nil
end

--- Calculates the dot product between this LuaVector object and another LuaVector object.
-- @param vector (table) The LuaVector object to calculate the dot product with.
-- @return (number) The dot product value.
function LuaVector:Dot(vector)
  return self[1] * vector[1] + self[2] * vector[2] + self[3] * vector[3]
end

--- Calculates the cross product between this LuaVector object and another LuaVector object.
-- @param vector (table) The LuaVector object to calculate the cross product with.
-- @return (table) The cross product as a new LuaVector object.
function LuaVector:Cross(vector)
  return LuaVector.new(
    self[2] * vector[3] - self[3] * vector[2],
    self[3] * vector[1] - self[1] * vector[3],
    self[1] * vector[2] - self[2] * vector[1]
  )
end

--- Calculates the distance between this LuaVector object and another LuaVector object.
-- @param vector (table) The LuaVector object to calculate the distance to.
-- @return (number) The distance between the LuaVectors.
function LuaVector:Distance(vector)
  local diff = self - vector
  return diff:Length()
end

--- Calculates the squared distance between this LuaVector object and another LuaVector object.
-- @param vector (table) The LuaVector object to calculate the squared distance to.
-- @return (number) The squared distance between the LuaVectors.
function LuaVector:DistToSqr(vector)
  local diff = self - vector
  return diff:LengthSqr()
end

--- Returns the negated version of this LuaVector object.
-- @return (table) The negated LuaVector object.
function LuaVector:GetNegated()
  self.cache = nil
  return LuaVector.new(-self[1], -self[2], -self[3])
end

--- Returns the normalized version of this LuaVector object.
-- @return (table) The normalized LuaVector object.
function LuaVector:GetNormalized()
  local length = self:Length()
  if length > 0 then
    return self / length
  else
    return LuaVector.new(0, 0, 0)
  end
end

--- Checks if two LuaVector objects are equal within a tolerance.
-- @param compare (table) The LuaVector object to compare to.
-- @param tolerance (number) The tolerance value. Default is FLT_EPSILON * 10.
-- @return (boolean) True if the LuaVector objects are equal within the tolerance, false otherwise.
function LuaVector:IsEqualTol(compare, tolerance)
  tolerance = tolerance or LuaVector.FLT_EPSILON

  return math.abs(self[1] - compare[1]) <= tolerance and
         math.abs(self[2] - compare[2]) <= tolerance and
         math.abs(self[3] - compare[3]) <= tolerance
end

--- Checks if the LuaVector object is zero (has coordinates (0, 0, 0)).
-- @return (boolean) True if the LuaVector object is zero, false otherwise.
function LuaVector:IsZero()
  return self[1] == 0 and self[2] == 0 and self[3] == 0
end

--- Sets all coordinates of the LuaVector object to zero (0, 0, 0).
function LuaVector:Zero()
  self[1] = 0
  self[2] = 0
  self[3] = 0
  self.cache = nil
end

--- Rotates the LuaVector object around the z-axis by the given rotation.
-- @param rotation (Angle) The rotation to apply.
function LuaVector:Rotate(rotation)
  local rad = math.rad(rotation.yaw)
  local cos = math.cos(rad)
  local sin = math.sin(rad)

  local x = self[1] * cos - self[2] * sin
  local y = self[1] * sin + self[2] * cos

  self[1] = x
  self[2] = y
  self.cache = nil
end

--- Sets the coordinates of the LuaVector object to match the given vector.
-- @param vector (table) The vector to set the coordinates from.
function LuaVector:Set(vector)
  self[1] = vector[1]
  self[2] = vector[2]
  self[3] = vector[3]
  self.cache = nil
end

--- Sets the coordinates of the LuaVector object to the given unpacked values.
-- @param x (number) The x-coordinate.
-- @param y (number) The y-coordinate.
-- @param z (number) The z-coordinate.
function LuaVector:SetUnpacked(x, y, z)
  self[1] = x
  self[2] = y
  self[3] = z
  self.cache = nil
end

--- Converts the LuaVector object to a table.
-- @return (table) The LuaVector as a table with keys 'x', 'y', and 'z'.
function LuaVector:ToTable()
  return self
end

--- Unpacks the coordinates of the LuaVector object.
-- @return (number, number, number) The unpacked coordinates (x, y, z).
function LuaVector:Unpack()
  return self[1], self[2], self[3]
end

--- Checks if the LuaVector object is within the axis-aligned box defined by the start and end points.
-- @param boxStart (table) The start point of the box as a LuaVector object.
-- @param boxEnd (table) The end point of the box as a LuaVector object.
-- @return (boolean) True if the LuaVector object is within the box, false otherwise.
function LuaVector:WithinAABox(boxStart, boxEnd)
  return self[1] >= boxStart[1] and self[2] >= boxStart[2] and self[3] >= boxStart[3] and
         self[1] <= boxEnd[1] and self[2] <= boxEnd[2] and self[3] <= boxEnd[3]
end

--- Converts the LuaVector object to a Garry's Mod Vector object. Cached for performance reasons until the vector is dirtied.
-- @return (Vector) The LuaVector as a Garry's Mod Vector object.
function LuaVector:ToGModVector()
  if not self.cache then
    self.cache = Vector(self[1], self[2], self[3])
  end
  return self.cache
end

--- Calculates the pitch, yaw, and roll angles from the LuaVector object.
-- @return (Angle) The calculated Angle object.
function LuaVector:Angle()
  local forward = self:GetNormalized()
  local angles = {}

  if forward[2] == 0 and forward[1] == 0 then
    angles.yaw = 0
    if forward[3] > 0 then
      angles.pitch = 270
    else
      angles.pitch = 90
    end
  else
    angles.yaw = (math.atan2(forward[2], forward[1]) * 180 / math.pi)
    if angles.yaw < 0 then
      angles.yaw = angles.yaw + 360
    end

    local tmp = math.sqrt(forward[1] * forward[1] + forward[2] * forward[2])
    angles.pitch = (math.atan2(-forward[3], tmp) * 180 / math.pi)
    if angles.pitch < 0 then
      angles.pitch = angles.pitch + 360
    end
  end

  angles.roll = 0

  return Angle(angles.pitch, angles.yaw, angles.roll)
end

--- Calculates the pitch, yaw, and roll angles from the LuaVector object with a specified up vector.
-- @param up (table) The up vector as a LuaVector object.
-- @return (Angle) The calculated Angle object.
function LuaVector:AngleEx(up)
  if (isvector(up)) then
    up = LuaVector(up)
  end

  local forward = self:GetNormalized()
  local left = up:Cross(forward):GetNormalized()
  local angles = {}

  local xyDist = math.sqrt(forward[1] * forward[1] + forward[2] * forward[2])

  if xyDist > 0.001 then
    angles.yaw = math.deg(math.atan2(forward[2], forward[1]))
    angles.pitch = math.deg(math.atan2(-forward[3], xyDist))

    local up_z = left[2] * forward[1] - left[1] * forward[2]
    angles.roll = math.deg(math.atan2(left[3], up_z))
  else
    angles.yaw = math.deg(math.atan2(-left[1], left[2]))
    angles.pitch = math.deg(math.atan2(-forward[3], xyDist))
    angles.roll = 0
  end

  return Angle(angles.pitch, angles.yaw, angles.roll)
end

--- Performs an OBB on OBB intersection test between two LuaVector objects.
-- @param box1Origin (table) The origin of the first box as a LuaVector object.
-- @param box1Angles (Angle) The angles of the first box as an Angle object.
-- @param box1Mins (table) The minimum bounds of the first box as a LuaVector object.
-- @param box1Maxs (table) The maximum bounds of the first box as a LuaVector object.
-- @param box2Origin (table) The origin of the second box as a LuaVector object.
-- @param box2Angles (Angle) The angles of the second box as an Angle object.
-- @param box2Mins (table) The minimum bounds of the second box as a LuaVector object.
-- @param box2Maxs (table) The maximum bounds of the second box as a LuaVector object.
-- @return (boolean) True if the boxes intersect, false otherwise.
function LuaVector.IsOBBIntersectingOBB(box1Origin, box1Angles, box1Mins, box1Maxs, box2Origin, box2Angles, box2Mins, box2Maxs)
  -- Convert the box angles to radians
  local box1RadAngles = box1Angles:ToRad()
  local box2RadAngles = box2Angles:ToRad()

  -- Compute the rotation matrices
  local box1RotMatrix = Matrix()
  box1RotMatrix:SetAngles(box1RadAngles)
  
  local box2RotMatrix = Matrix()
  box2RotMatrix:SetAngles(box2RadAngles)

  -- Compute the center points of the boxes
  local box1Center = (box1Mins + box1Maxs) * 0.5
  local box2Center = (box2Mins + box2Maxs) * 0.5

  -- Compute the extents of the boxes
  local box1Extents = (box1Maxs - box1Mins) * 0.5
  local box2Extents = (box2Maxs - box2Mins) * 0.5

  -- Compute the axes for the boxes
  local box1AxisX = box1RotMatrix:GetForward()
  local box1AxisY = box1RotMatrix:GetRight()
  local box1AxisZ = box1RotMatrix:GetUp()

  local box2AxisX = box2RotMatrix:GetForward()
  local box2AxisY = box2RotMatrix:GetRight()
  local box2AxisZ = box2RotMatrix:GetUp()

  -- Check for separation along each axis
  for _, axis in ipairs({ box1AxisX, box1AxisY, box1AxisZ, box2AxisX, box2AxisY, box2AxisZ }) do
    local box1Projection = box1Extents:Dot(axis)
    local box2Projection = box2Extents:Dot(axis)

    local dist = math.abs((box2Origin - box1Origin):Dot(axis))

    if dist > box1Projection + box2Projection then
      -- Boxes are separated along this axis, they don't intersect
      return false
    end
  end

  -- Check for separation along the cross products of box axes
  for _, axis1 in ipairs({ box1AxisX, box1AxisY, box1AxisZ }) do
    for _, axis2 in ipairs({ box2AxisX, box2AxisY, box2AxisZ }) do
      local axis = axis1:Cross(axis2)
      local box1Projection = box1Extents:Dot(axis)
      local box2Projection = box2Extents:Dot(axis)

      local dist = math.abs((box2Origin - box1Origin):Dot(axis))

      if dist > box1Projection + box2Projection then
        -- Boxes are separated along this axis, they don't intersect
        return false
      end
    end
  end

  -- Boxes overlap along all axes, they intersect
  return true
end

--- Calculates the intersection point and parameter along a ray with an axis-aligned box (OBB).
-- @param rayStart (table) The start point of the ray as a LuaVector object.
-- @param rayDelta (table) The direction and length of the ray as a LuaVector object.
-- @param boxOrigin (table) The origin of the box as a LuaVector object.
-- @param boxAngles (Angle) The angles of the box as an Angle object.
-- @param boxMins (table) The minimum bounds of the box as a LuaVector object.
-- @param boxMaxs (table) The maximum bounds of the box as a LuaVector object.
-- @return (table, number) The intersection point as a LuaVector object and the parameter along the ray.
function LuaVector:IntersectRayWithOBB(rayStart, rayDelta, boxOrigin, boxAngles, boxMins, boxMaxs)
  local rayStartRelativeToBox = rayStart - boxOrigin
  rayStartRelativeToBox:Rotate(boxAngles:Inverse())

  local tStart = 0
  local tEnd = 1

  for i = 1, 3 do
    if math.abs(rayDelta[i]) < LuaVector.FLT_EPSILON then
      if rayStartRelativeToBox[i] < boxMins[i] or rayStartRelativeToBox[i] > boxMaxs[i] then
        return nil
      end
    else
      local t1 = (boxMins[i] - rayStartRelativeToBox[i]) / rayDelta[i]
      local t2 = (boxMaxs[i] - rayStartRelativeToBox[i]) / rayDelta[i]

      if t1 > t2 then
        t1, t2 = t2, t1
      end

      if t1 > tStart then
        tStart = t1
      end
      if t2 < tEnd then
        tEnd = t2
      end

      if tStart > tEnd then
        return nil
      end
    end
  end

  local hitPos = rayStart + rayDelta * tStart
  hitPos:Rotate(boxAngles)

  return hitPos, tStart
end

--- Calculates the intersection point between a ray and a plane defined by a position and normal.
-- @param rayOrigin (table) The origin of the ray as a LuaVector object.
-- @param rayDirection (table) The direction of the ray as a LuaVector object.
-- @param planePosition (table) The position of the plane as a LuaVector object.
-- @param planeNormal (table) The normal vector of the plane as a LuaVector object.
-- @return (table) The intersection point as a LuaVector object, or nil if there is no intersection.
function LuaVector:IntersectRayWithPlane(rayOrigin, rayDirection, planePosition, planeNormal)
  local denom = self:Dot(planeNormal)
  if math.abs(denom) < LuaVector.FLT_EPSILON then
    return nil
  end

  local t = (planePosition - rayOrigin):Dot(planeNormal) / denom

  if t >= 0 then
    return rayOrigin + rayDirection * t
  else
    return nil
  end
end

--- Calculates the distance between a line defined by two points and a point in 3D space.
-- @param lineStart (table) The start point of the line as a LuaVector object.
-- @param lineEnd (table) The end point of the line as a LuaVector object.
-- @param pointPos (table) The position of the point as a LuaVector object.
-- @return (number) The distance between the line and the point.
function LuaVector:DistanceToLine(lineStart, lineEnd, pointPos)
  local lineDir = (lineEnd - lineStart):GetNormalized()
  local lineToPoint = pointPos - lineStart

  local dot = lineToPoint:Dot(lineDir)
  local closestPoint = lineStart + lineDir * dot

  return (pointPos - closestPoint):Length()
end

setmetatable(LuaVector, { 
  __call = function (self, ...)
    return LuaVector.new(...) 
  end
})

setmetatable(LuaVector, { 
  __call = function (self, ...)
    return LuaVector.new(...) 
  end
})

local function testHelper(x, y, z, func, ...)
  local LuaVec = LuaVector(x, y, z)
  local GmodVec = Vector(x, y, z)

  local r1 = { LuaVector[func](LuaVec, ...) }
  local r2 = { GmodVec[func](GmodVec, ...) }

  assert(#r1 == #r2, func .. ": Expected " .. #r2 .. " return values, got " .. #r1)

  for k, v in pairs(r1) do
    local v2 = r2[k]

    local t = type(v)
    local succ
    if (t == "number") then
      succ = (v - v2) < LuaVector.FLT_EPSILON
    elseif (t == "boolean") then
      succ = v == v2
    elseif (t == "Angle") then
      v = LuaVector.NormalizeAngle(v)
      v2 = LuaVector.NormalizeAngle(v2)
      succ = v == v2
    elseif v.__type == "LuaVector" then
      succ = v:IsEqualTol(v2)
    else
      succ = v == v2
    end
    assert(succ, func .. "[" .. k .. "]: Expected (" .. t .. ") " .. tostring(v) .. ", got (" .. type(v2) .. ") " .. tostring(v2))
  end

  return r1, r2
end

function LuaVector.NormalizeAngle(ang)
  return Angle(math.Truncate(math.NormalizeAngle(ang.pitch), LuaVector.NORMALIZE_ANGLE_DP), math.Truncate(math.NormalizeAngle(ang.yaw), LuaVector.NORMALIZE_ANGLE_DP), math.Truncate(math.NormalizeAngle(ang.roll), LuaVector.NORMALIZE_ANGLE_DP))
end

function LuaVector.RunTest(x, y, z)
  testHelper(x, y, z, "Add", Vector(x^2, y^2, z^2))
  testHelper(x, y, z, "Sub", Vector(x^2, y^2, z^2))
  testHelper(x, y, z, "Mul", Vector(x^2, y^2, z^2))
  testHelper(x, y, z, "Div", x^2)
  testHelper(x, y, z, "Length")
  testHelper(x, y, z, "LengthSqr")
  testHelper(x, y, z, "Normalize")
  testHelper(x, y, z, "Dot", Vector(x^2, y^2, z^2))
  testHelper(x, y, z, "Cross", Vector(x^2, y^2, z^2))
  testHelper(x, y, z, "Distance", Vector(x^2, y^2, z^2))
  testHelper(x, y, z, "DistToSqr", Vector(x^2, y^2, z^2))
  testHelper(x, y, z, "GetNegated")

  testHelper(x, y, z, "IsZero")
  testHelper(0, 0, 0, "IsZero")

  testHelper(x, y, z, "Zero")
  testHelper(x, y, z, "Rotate", Angle(x^2, y^2, z^2))
  testHelper(x, y, z, "Set", Vector(x^2, y^2, z^2))
  testHelper(x, y, z, "SetUnpacked", x^2, y^2, z^2)
  testHelper(x, y, z, "Unpack")

  -- Inside case
  local boxStart1, boxEnd1 = Vector(x-1, y-1, z-1), Vector(x+1, y+1, z+1)
  testHelper(x, y, z, "WithinAABox", boxStart1, boxEnd1)
  
  -- Outside case
  local boxStart2, boxEnd2 = Vector(x+1, y+1, z+1), Vector(x+2, y+2, z+2)
  testHelper(x, y, z, "WithinAABox", boxStart2, boxEnd2)

  -- Boundary case 1
  local boxStart2, boxEnd2 = Vector(x, y, z), Vector(x+2, y+2, z+2)
  testHelper(x, y, z, "WithinAABox", boxStart2, boxEnd2)

  -- Boundary case 2
  local boxStart2, boxEnd2 = Vector(x-1, y-1, z-1), Vector(x, y, z)
  testHelper(x, y, z, "WithinAABox", boxStart2, boxEnd2)

  testHelper(x, y, z, "Angle")
  testHelper(x, y, z, "AngleEx", Vector(x^2, y^2, z^2):GetNormalized())
  

  -- IsEqualTol true case
  local lv = LuaVector(x, y, z)
  local gv = Vector(x, y, z)
  assert(lv:IsEqualTol(gv), "IsEqualTol: Expected " .. tostring(gv) .. ", got " .. tostring(lv))

  -- False case
  gv = Vector(x^2, y^2, z^2)
  assert(not lv:IsEqualTol(gv), "IsEqualTol: Expected not to be equal = " .. tostring(gv) .. ", got " .. tostring(lv))

  -- Make sure tables match
  gv = Vector(x, y, z)
  lv = LuaVector(x, y, z)
  local gv_t = gv:ToTable()
  local lv_t = lv:ToTable()

  for k, v in pairs(gv_t) do
    assert(v == lv_t[k], "ToTable: Expected " .. tostring(v) .. ", got " .. tostring(lv_t[k]))
  end

  return true
end

return LuaVector
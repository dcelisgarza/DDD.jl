"""
```
```
"""
function makeSegment(
    type::AbstractDlnSeg,
    slipPlane::Vector{T},
    bVec::Vector{T},
) where {T<:Real}
    if typeof(type) == segEdge
        seg = cross(slipPlane, bVec)
    elseif typeof(type) == segEdgeN
        seg = slipPlane
    elseif typeof(type) == segScrew
        seg = bVec
    elseif typeof(type) == segMixed
        error("makeSegment: mixed dislocation segment undefined")
    elseif typeof(type) == segNone
        return [0.0; 0.0; 0.0]
    end
    return seg ./= norm(seg)
end
"""
```
```
Idealised dislocation loop.
"""
struct DislocationLoop{
    T1<:loopSides,
    T2<:Integer,
    T3<:Union{
        T8 where {T8<:AbstractDlnSeg},
        AbstractArray{<:AbstractDlnSeg,N} where {N},
    },
    T4<:Union{T9 where {T9<:Real},AbstractArray{<:Real,N} where {N}},
    T5<:Union{Integer,AbstractArray{<:Integer,N} where {N}},
    T6<:AbstractArray{<:Integer,N} where {N},
    T7<:AbstractArray{<:Real,N} where {N},
    T8<:Vector{<:nodeType},
}
    numSides::T1
    nodeSide::T2
    segType::T3
    segLen::T4
    slipSystem::T5
    links::T6
    slipPlane::T7
    bVec::T7
    coord::T7
    label::T8
    numLoops::T2
    function DislocationLoop(
        numSides,
        nodeSide,
        segType,
        segLen,
        slipSystem,
        _slipPlane,
        _bVec,
        label,
        numLoops,
    )
        nodeTotal = numSides * nodeSide
        numSegType = length(segType)
        @assert length(label) == nodeTotal
        @assert numSegType == length(segLen) == numSides / 2 ==
                size(_bVec, 1) == size(_slipPlane, 1) == length(slipSystem)
        links = zeros(Integer, nodeTotal, 2)
        coord = zeros(nodeTotal, 3)
        seg = zeros(numSegType, 3)
        slipPlane = zeros(0, 3)
        bVec = zeros(0, 3)

        for i = 1:numSegType
            seg[i, :] = makeSegment(segType[i], _slipPlane[i, :], _bVec[i, :]) *
                        segLen[i]
        end
        if numSides == 4
            # First node.
            coord[1, :] -= seg[1, :] + seg[2, :]
            slipPlane = [slipPlane; _slipPlane[2, :]']
            bVec = [bVec; _bVec[2, :]']
            for k = 1:nodeSide
                coord[1+k, :] = coord[k, :] + seg[1, :]
                slipPlane = [slipPlane; _slipPlane[1, :]']
                bVec = [bVec; _bVec[1, :]']
            end
            for k = 1:nodeSide
                coord[1 + k + nodeSide, :] = coord[k+nodeSide, :] + seg[2, :]
                slipPlane = [slipPlane; _slipPlane[2, :]']
                bVec = [bVec; _bVec[2, :]']
            end
            for k = 1:nodeSide
                coord[1 + k + 2 * nodeSide, :] = coord[k+2*nodeSide, :] -
                                                 seg[1, :]
                slipPlane = [slipPlane; _slipPlane[1, :]']
                bVec = [bVec; _bVec[1, :]']
            end
            if nodeSide > 1
                # Last node
                coord[2+3*nodeSide, :] = coord[1+3*nodeSide, :] - seg[2, :]
                slipPlane = [slipPlane; _slipPlane[2, :]']
                bVec = [bVec; _bVec[2, :]']
            end
        elseif numSides == 6
            # First node.
            coord[1, :] -= seg[1, :] + seg[2, :] + seg[3, :]
            slipPlane = [slipPlane; _slipPlane[3, :]']
            bVec = [bVec; _bVec[3, :]']
            for k = 1:nodeSide
                coord[1+k, :] = coord[k, :] + seg[1, :]
                slipPlane = [slipPlane; _slipPlane[1, :]']
                bVec = [bVec; _bVec[1, :]']
            end
            for k = 1:nodeSide
                coord[1 + k + nodeSide, :] = coord[k+nodeSide, :] + seg[2, :]
                slipPlane = [slipPlane; _slipPlane[2, :]']
                bVec = [bVec; _bVec[2, :]']
            end
            for k = 1:nodeSide
                coord[1 + k + 2 * nodeSide, :] = coord[k+2*nodeSide, :] +
                                                 seg[3, :]
                slipPlane = [slipPlane; _slipPlane[3, :]']
                bVec = [bVec; _bVec[3, :]']
            end
            for k = 1:nodeSide
                coord[1 + k + 3 * nodeSide, :] = coord[k+3*nodeSide, :] -
                                                 seg[1, :]
                slipPlane = [slipPlane; _slipPlane[1, :]']
                bVec = [bVec; _bVec[1, :]']
            end
            for k = 1:nodeSide
                coord[1 + k + 4 * nodeSide, :] = coord[k+4*nodeSide, :] -
                                                 seg[2, :]
                slipPlane = [slipPlane; _slipPlane[2, :]']
                bVec = [bVec; _bVec[2, :]']
            end
            if nodeSide > 1
                for k = 1:nodeSide-1
                    coord[1 + k + 5 * nodeSide, :] = coord[k+5*nodeSide, :] -
                                                     seg[3, :]
                    slipPlane = [slipPlane; _slipPlane[3, :]']
                    bVec = [bVec; _bVec[3, :]']
                end
            end
        else
            error("more sides for a source loop are undefined")
        end
        # Links
        for j = 1:nodeTotal-1
            links[j, :] = [j; 1 + j]
        end
        links[nodeTotal, :] = [nodeTotal; 1]

        new{
            typeof(numSides),
            typeof(nodeSide),
            typeof(segType),
            typeof(segLen),
            typeof(slipSystem),
            typeof(links),
            typeof(slipPlane),
            typeof(label),
        }(
            numSides,
            nodeSide,
            segType,
            segLen,
            slipSystem,
            links,
            slipPlane,
            bVec,
            coord,
            label,
            numLoops,
        )
    end
end

function zero(::Type{DislocationLoop})
    DislocationLoop(
        loopSides(4),
        1,
        [segNone(), segNone()],
        [0.0, 0.0],
        zeros(Integer, 2),
        zeros(2, 3),
        zeros(2, 3),
        zeros(nodeType, 4),
        0,
    )
end
mutable struct DislocationNetwork{
    T1<:AbstractArray{<:Integer,N} where {N},
    T2<:AbstractArray{<:Real,N} where {N},
    T3<:Vector{nodeType},
    T4<:Integer,
}
    links::T1 # Links.
    slipPlane::T2 # Slip planes.
    bVec::T2 # Burgers vectors.
    coord::T2 # Node coordinates.
    label::T3 # Node labels.
    numNode::T4 # Number of dislocations.
    numSeg::T4 # Number of segments.
    function DislocationNetwork(
        links,
        slipPlane,
        bVec,
        coord,
        label,
        numNode = 0,
        numSeg = 0,
    )
        @assert size(links, 2) == 2
        @assert size(bVec, 2) == size(slipPlane, 2) == size(coord, 2) == 3
        @assert size(links, 1) == size(bVec, 1) == size(slipPlane, 1)
        @assert size(coord, 1) == size(label, 1)
        new{typeof(links),typeof(bVec),typeof(label),typeof(numNode)}(
            links,
            slipPlane,
            bVec,
            coord,
            label,
            numNode,
            numSeg,
        )
    end # Constructor
end # DislocationNetwork

function vcat(network::DislocationNetwork, n::Integer)
    network.links = [network.links; zeros(Integer, n, 2)]
    network.slipPlane = [network.slipPlane; zeros(n, 3)]
    network.bVec = [network.bVec; zeros(n, 3)]
    network.coord = [network.coord; zeros(n, 3)]
    network.label = [network.label; zeros(nodeType, n)]
    return network
end






struct DislocationP{T1<:Real,T2<:Integer,T3<:Bool,T4<:AbstractMobility}
    # Size.
    coreRad::T1 # Core radius.
    coreRadMag::T1 # Magnitude of core Radius.
    # Connectivity.
    minSegLen::T1 # Minimum line length.
    maxSegLen::T1 # Maximum line length.
    minArea::T1 # Minimum area for remeshing.
    maxArea::T1 # Maximum area for remeshing.
    maxConnect::T2 # Maximum number of connections to a node.
    remesh::T3 # Flag for remeshing.
    collision::T3 # Flag for collision handling.
    separation::T3 # Flag for separation handling.
    virtualRemesh::T3 # Flag for virtual remeshing.
    # Mobility.
    edgeDrag::T1 # Drag coefficient edge dislocation.
    screwDrag::T1 # Drag coefficient screw dislocation.
    climbDrag::T1 # Drag coefficient climb.
    lineDrag::T1 # Drag coefficient line.
    mobility::T4 # Mobility law.
    # Fool-proof constructor.
    function DislocationP(
        coreRad,
        coreRadMag,
        minSegLen,
        maxSegLen,
        minArea,
        maxArea,
        maxConnect,
        remesh,
        collision,
        separation,
        virtualRemesh,
        edgeDrag,
        screwDrag,
        climbDrag,
        lineDrag,
        mobility,
    )
        coreRad == minSegLen == maxSegLen == 0 ? nothing :
        @assert coreRad < minSegLen < maxSegLen
        minArea == maxArea == 0 ? nothing : @assert minArea < maxArea
        new{
            typeof(coreRad),
            typeof(maxConnect),
            typeof(remesh),
            typeof(mobility),
        }(
            coreRad,
            coreRadMag,
            minSegLen,
            maxSegLen,
            minArea,
            maxArea,
            maxConnect,
            remesh,
            collision,
            separation,
            virtualRemesh,
            edgeDrag,
            screwDrag,
            climbDrag,
            lineDrag,
            mobility,
        )
    end # constructor
end # DislocationP

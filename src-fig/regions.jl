# Just remembering the masks for each region
# Clockwise, without end point
regions = (;
    arrest = [
        Point2f(-1750, -30),
        Point2f(-2000, -3),
        Point2f(-1700, -3),
        Point2f(-1050, -30),
        Point2f(-1050, -55),
    ],
    total = [
        Point2f(-2500, -3), 
        Point2f(2500, -3), 
        Point2f(2500, -100), 
        Point2f(-2500, -100),
    ]
)
region_names = (; 
    arrest="Arrest region",
    total="Whole front"
)
#=
    top = [
        Point2f(1250, -1),
        Point2f(1250, -80),
        Point2f(250, -52),
        Point2f(-1000, -25),
        Point2f(-1650, -1),
    ],
    entrainment = [
        Point2f(-900, -25),
        Point2f(250, -58),
        Point2f(1250, -83),
        Point2f(1250, -93),
        Point2f(750, -83),
        Point2f(-900, -38),
    ],
    underneath = [
        Point2f(750, -90),
        Point2f(-900, -90),
        Point2f(-900, -45),
        Point2f(750, -85),
    ]

    top="Top",
    entrainment="Entrainment",
    underneath="Underneath",
=#
linecheck(x, z, p₁, p₂) = linecheck(x, z, p₁[1], p₁[2], p₂[1], p₂[2])
linecheck(x, z, x₁, z₁, x₂, z₂) = (x - x₁) * (z₂ - z₁) - (z - z₁) * (x₂ - x₁) < 0
maskfromlines(x, z, points) = mapreduce(&, points, circshift(points, 1)) do p₁, p₂
    linecheck(x, z, p₁, p₂)
end

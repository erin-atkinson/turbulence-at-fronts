function subfig_label!(scene, label; fontsize=14)
    gl = GridLayout(scene, 
        tellwidth = false, 
        tellheight = false, 
        halign = :left, 
        valign = :top,
    )
    Box(gl[1, 1], color = :white, strokecolor = :black, strokewidth = 0)
    Label(gl[1, 1], label;
        fontsize,
        padding = (1, 1, 2, 2),
    )
end

# Pass an integer to get that letter 
function subfig_label!(scene, label::Integer; fontsize=14)
    char = 'a' + (label - 1)
    subfig_label!(scene, "$char)"; fontsize)
end

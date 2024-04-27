# 该文件旨在设定工作路径
wkdir="/Users/kangjiayin/Desktop/julia/DSE-0419-ps/"
cd(wkdir)

function modify_file_line(filepath::AbstractString, line_num::Integer, new_content::AbstractString)
    # 读取文件
    f = open(filepath, "r")
    lines = readlines(f)
    close(f)
    
    # 修改指定行的内容
    lines[line_num] = new_content
    
    # 写回到文件中
    f = open(filepath, "w")
    foreach(line -> write(f, line, "\n"), lines)
    close(f)
end

for epsilon = 0.2:0.1:0.2
    for omega2num = 1.1:0.5:1.1
    modify_file_line("./config.toml",46,"epsilon = $epsilon")
    modify_file_line("./config.toml",47,"omega2 = $omega2num")
    include(joinpath(pwd(),"src/main.jl"))
    end 
end

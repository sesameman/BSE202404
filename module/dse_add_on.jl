module dse_add_on

export print_seconds

function print_seconds(seconds)
    if seconds < 60
        seconds = round(seconds, sigdigits=4)
        println("Total time: $seconds seconds")
    elseif seconds < 3600
        seconds /= 60
        seconds = round(seconds, sigdigits=4)
        println("Total time: $seconds minutes")
    elseif seconds < 86400
        seconds /= 3600
        seconds = round(seconds, sigdigits=4)
        println("Total time: $seconds hours")
    else
        seconds /= 86400
        seconds = round(seconds, sigdigits=4)
        println("Total time: $seconds days")
    end
end

end
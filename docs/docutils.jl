using Pkg
using Revise
using Printf
using DelimitedFiles


function extract_function_names(file_path::AbstractString)
    # Read the content of the file
    content = read(file_path, String)

    # Extract function names using a regular expression
    function_names_itterator =eachmatch(r"(function\s+(([a-zA-z_]*)\w*(\(.*\))))", content);
    function_names =  [aa.captures[2] for aa in function_names_itterator];
    return function_names
end

function extract_docstrings(file_path::AbstractString)
    # Read the content of the file
    content = read(file_path, String)
    function_docstr = Dict{String,String}()
    # Extract function names using a regular expression
    reg=r"\"\"\"\n\s*(.*)\s*DOCSTRING\n([\s\S]*?)\"\"\"";
    matches = [ea.captures for ea in eachmatch(reg,content)];
    matches = string.(permutedims(hcat(matches...)))
    return matches
end

function extract_regex(file_path::AbstractString,reg)
    # Read the content of the file
    content = read(file_path, String)
    function_docstr = Dict{String,String}()
    # Extract function names using a regular expression
    function_names_itterator =eachmatch(reg, content);
    function_names_obj =  [aa for aa in function_names_itterator];
    return function_names_obj
end

function library_docstrings(file_path::AbstractString)
    # Read the content of the file
    content = read(file_path, String)
    reg = r"\"\"\"\n\s*(.*)\s*DOCSTRING\n([\s\S]*?)ELEM TYPE:\s*([A-Z]*)\nEQUATIONS:\n([\s\S]*?)INPUTS\n([\s\S]*?)\"\"\"";

    # captures
    matches = [ea.captures for ea in eachmatch(reg,content)];
    matches = string.(permutedims(hcat(matches...)))

    return matches
end


# ## generates text files for component doc strings
# function generate_libtxt(filename::String,docstrings::Matrix{String})
#     txtpath = "docs/pre/" * filename * ".txt";
#     # Tpes of components
#     ELEM_TYPES   = docstrings[:,3]
#     UNIQUE_ELEM  = unique(ELEM_TYPES)
#     ORDER        = sortperm(UNIQUE_ELEM);
#     HEADERS = UNIQUE_ELEM[ORDER];
#     html = "<h1>"*filename*"<h1>\n"

#     prep(x) = replace(strip(x),r"\n\s*"=>"\n");

#     for hd in HEADERS
#         html *= "\n\n<div><section><h1>"*hd*"</h1>"
#         fcn_idx = findall(x -> x==hd, docstrings[:,3])
#         for idx in fcn_idx
#             row = docstrings[idx,:];
#             html *= "\n\t<h4>"*row[1]*"</h4>\n"
#             html *= "\t\t<ul>\n"
#             html *= "\t\t<li>\n\t\t\t<h5>DESCRIPTION:</h5>\n"
#             html *= "\t\t\t<p>\n\t\t\t\t" * replace(prep(row[2]),"\n" => "<br>\n\t\t\t")* "\n\t\t\t</p>\n\t\t</li>\n"
#             html *= "\t\t<li>\n\t\t\t<h5>INPUTS:</h5>\n"
#             html *= "\t\t\t<ul>\n"
#             html *= "\t\t\t\t<li>" * replace(row[5],"\n" => "</li>\n\t\t\t\t<li>") * "\t\t\t</li>\n\t\t\t</ul>\n\t\t</li>\n"
#             html *= "\t\t<li>\n\t\t\t<h5>EQUATIONS:</h5>\n"
#             html *= "\t\t\t\t<p>\n\t\t\t\t\\[ " * replace(replace(row[4],r"#.*\n" => "\n\t\t\t\t"),"\n" => "\\]\n\t\t\t\t\\[ " ,"~" => "=")  * "\\]\n\t\t</li>\n"
#             html *= "\t\t\t\t</p>\n\t\t\t</ul>\n"
#         end
#         html *= "\n</section></div>\n"
#     end
#     touch(txtpath)
#     write(txtpath,html)
#     # writedlm(txtpath,html,quotes=false)
#     # open(txtpath,"w") do file
#     #     write(file,html)
#     # end
# end

function generate_libtxt(filename::String,file_path::AbstractString)
    docstrings = library_docstrings(file_path);
    fcn_names = extract_function_names(file_path);
    fcn_names = fcn_names[sortperm(fcn_names)];
    txtpath = "docs/pre/" * filename * ".txt";
    # Tpes of components
    ELEM_TYPES   = docstrings[:,3]
    UNIQUE_ELEM  = unique(ELEM_TYPES)
    ORDER        = sortperm(UNIQUE_ELEM);
    HEADERS = UNIQUE_ELEM[ORDER];
    html = "<h1>"*filename*"<h1>\n"
    html *= "\n\n<div><section><h3> All Functions </h3>\n"
    html *= "\t\t<ul>\n"

    for fn in fcn_names
        html *= "\t\t\t<li>" * fn * "</li>\n"
    end
    html *= "\t\t</ul>\n"


    prep(x) = replace(strip(x),r"\n\s*"=>"\n");

    for hd in HEADERS
        html *= "\n\n<div><section><h3>"*hd*"S"*"</h3>"
        fcn_idx = findall(x -> x==hd, docstrings[:,3])
        for idx in fcn_idx
            row = docstrings[idx,:];
            html *= "\n\t<h4>"*row[1]*"</h4>\n"
            html *= "\t\t<ul>\n"
            html *= "\t\t<li>\n\t\t\t<h5>DESCRIPTION:</h5>\n"
            html *= "\t\t\t<p>\n\t\t\t\t" * replace(prep(row[2]),"\n" => "<br>\n\t\t\t")* "\n\t\t\t</p>\n\t\t</li>\n"
            html *= "\t\t<li>\n\t\t\t<h5>INPUTS:</h5>\n"
            html *= "\t\t\t<ul>\n"
            html *= "\t\t\t\t<li>" * replace(row[5],"\n" => "</li>\n\t\t\t\t<li>") * "\t\t\t</li>\n\t\t\t</ul>\n\t\t</li>\n"
            html *= "\t\t<li>\n\t\t\t<h5>EQUATIONS:</h5>\n"
            html *= "\t\t\t\t<p>\n\t\t\t\t\\[ " * replace(replace(row[4],r"#.*\n" => "\n\t\t\t\t"),"\n" => "\\]\n\t\t\t\t\\[ " ,"~" => "=")  * "\\]\n\t\t</li>\n"
            html *= "\t\t\t\t</p>\n\t\t\t</ul>\n"
        end
        html *= "\n</section></div>\n"
    end
    touch(txtpath)
    write(txtpath,html)
end

function generate_utiltxt(filename::String,file_path::AbstractString)
    docstrings = extract_docstrings(file_path);
    fcn_names = extract_function_names(file_path);
    fcn_names = fcn_names[sortperm(fcn_names)];
    docstrings = docstrings[sortperm(fcn_names),:];
    txtpath = "docs/pre/" * filename * ".txt";
    # Tpes of components
    html = "<h1>"*filename*"<h1>\n"
    html *= "\n\n<div><section><h1> All Functions </h1>\n"
    html *= "\t\t<ul>\n"

    for fn in fcn_names
        html *= "\t\t\t<li>" * fn * "</li>\n"
    end
    html *= "\t\t</ul>\n"


    prep(x) = replace(strip(x),r"\n\s*"=>"\n");
    HEADERS = ["Functions"]
    for hd in HEADERS
        html *= "\n\n<div><section><h3>"*hd*"</h3>"
        for idx = 1:length(fcn_names)
            row = docstrings[idx,:];
            html *= "\n\t<h4>"*row[1]*"</h4>\n"
            html *= "\t\t<ul>\n"
            html *= "\t\t<li>\n\t\t\t<h5>DESCRIPTION:</h5>\n"
            html *= "\t\t\t<p>\n\t\t\t\t" * replace(prep(row[2]),"\n" => "<br>\n\t\t\t")* "\n\t\t\t</p>\n\t\t</li>\n"
            html *= "\t\t\t</ul>\n"
        end
        html *= "\n</section></div>\n"
    end
    touch(txtpath)
    write(txtpath,html)
end
"""     
    buildhtml(filref)


    DOCSTRING
    fileref is a string and fileref.txt is a file in pre/
    Builds documentation html page from txt files in pre/
"""
function buildhtml(fileref)
    @printf "Building html file for : %s \n" fileref
    
    # @printf "Active directory : %s \n" splitdir(pwd())[end]

    prefixpath  = "docs/tmpl/prefix.txt";
    navpath     = "docs/tmpl/navtxt.txt";
    suffixpath  = "docs/tmpl/suffix.txt";

    fpath = "docs/pre/" * fileref * ".txt";
    htmlpath = "docs/" *fileref * ".html";

    allfiles = readdir("docs/pre/");
    @assert fileref*".txt" ∈ readdir("docs/pre/") "File not in path"

    prefix=readlines(prefixpath);
    nav = readlines(navpath);
    body = readlines(fpath);
    head = body[1];
    main = body[2:end]
    suff = readlines(suffixpath);

    outtext = vcat(prefix,head,nav,main,suff)
    touch(htmlpath)
    open(htmlpath,"w") do file
        writedlm(file,outtext,quotes=false)
    end
    return outtext
end

function extract_docstrings2(file_path)
    docstrings = Dict{String, Dict{String, Vector{String}}}()

    current_function = ""
    current_elem_type = ""
    current_section = ""

    open(file_path) do file
        for line in eachline(file)
            if occursin("""\"\"\"""", line)
                current_section = ""
                continue
            end

            if occursin("DOCSTRING", line) || occursin("ELEM TYPE", line) || occursin("EQUATIONS", line) || occursin("INPUTS", line)
                current_section = strip(split(line, ":")[1])
                if current_section == "ELEM TYPE"
                    current_elem_type = strip(split(line, ":")[2])
                    get!(docstrings, current_elem_type, Dict{String, Vector{String}}())
                end
                continue
            end

            if current_section != ""
                push!(get!(get!(docstrings, current_elem_type), current_function, Vector{String}()), line)
            elseif length(strip(line)) > 0
                current_function = strip(line)
            end
        end
    end

    return docstrings
end

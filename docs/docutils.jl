using Pkg

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
    function_names_itterator =eachmatch(r"\"\"\"\n\s*(.*)\s*\n([\s\S]*?)\"\"\"[\n|\r](function\s+(([a-zA-z_]*)\w*(\(.*\))))", content);
    function_names_obj =  [aa for aa in function_names_itterator];
    function_docstr =  Dict(String(fn.captures[4]) => String(fn.captures[2]) for fn in function_names_obj)
    return function_docstr
end
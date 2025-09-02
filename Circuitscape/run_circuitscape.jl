using  Circuitscape
folder_path="//home/xiongh/Biomod2/new/result_1/IBR/"
files = readdir(folder_path,join=true) |> x -> filter(endswith(".ini"),x)
for file in files
	println(file)
	compute(file)
end


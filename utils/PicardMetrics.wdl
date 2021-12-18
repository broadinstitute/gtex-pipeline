version 1.0

task LoadPicardMetrics{
	input {
		File picard_metrics
	}
	Int disk_space = 10 + ceil(size(picard_metrics))

	command <<<
		python <<-EOF 
		import json

		with open("~{picard_metrics}","r") as metrics:
			firstline=True
			line_count=0
			data={}
			for line in metrics:
				line=line.rstrip("\n")
				if line.startswith("#") or line.strip()=="":
					continue
				if firstline:
					headers=[*map(lambda x: x.strip(),line.split("\t"))]
					print("found headers: " + headers)
					for header in headers:
						data[header]=[]
					firstline=False
				else:
					for header,value in zip(headers,line.split("\t")):
						data[header].append(value)
						line_count+=1
					
		with open("dump.json","w") as out:
			json.dump(data,out)

		EOF
	>>>
	output {
		Map[String,Array[String]] values_map= read_json("dump.json")
		File json_dump = "dump.json"
	}

	runtime {
		docker: "python:latest"
		memory: "2GB"
		disks: "local-disk ~{disk_space} HDD"
	}
}

workflow LoadPicardMetricsWF{
	input {
		File picard_metrics
		String requested_metric
		Int requested_row
	}
	call LoadPicardMetrics{input: picard_metrics=picard_metrics}

	output {
		String requested_value=LoadPicardMetrics.values_map[requested_metric][requested_row]
		Map[String,Array[String]] values_map= LoadPicardMetrics.values_map
		File json_dump = LoadPicardMetrics.json_dump
	}
}
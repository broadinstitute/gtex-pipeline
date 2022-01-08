version 1.0 


workflow GetString{
	input {
		
		String string
	}
	#[["membership:sample_set_id","sample"]]+

	output {
		File out=string
	}	
}

version 1.0 


workflow CreateParticipantSampleSets{
	input {
		Array[String] samples
		Array[String] participants
	}
	#[["membership:sample_set_id","sample"]]+
	Array[Array[String]] membership_array =  flatten([[["membership:sample_set_id","sample"]] , transpose([participants,samples])])
	output {
		File memberships=write_tsv(membership_array)
	}	
}

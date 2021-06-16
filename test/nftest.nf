params.reads = "fake{1,2}.fastq"
params.R1_barcode = 'true'

test_ch = Channel.fromPath( params.reads )

reads =  test_ch.toSortedList()

println reads

if ( params.R1_barcode == 'true' ) {

myFile = new File(reads[0])
String line = myFile.readLines().get(1)
println line.length()
params.bc = 5
params.umi = 3
params.R1_barcode = 'true'
fulllen = params.bc + params.umi
difference = line.length() - fulllen
println line.length() - fulllen

}

params.reads = "fake{1,2}.fastq"
params.R1_barcode = 'true'
params.bc = 5
params.umi = 5

test_ch = Channel.fromPath( params.reads )

reads =  test_ch.toSortedList().get()

if ( params.R1_barcode == 'true' ) {

myFile = new File("${reads[0]}")
String line = myFile.readLines().get(1)
fulllen = params.bc + params.umi
difference = line.length() - fulllen

}

else if ( params.R1_barcode != 'true' ) {

myFile = new File("${reads[1]}")
String line = myFile.readLines().get(1)
fulllen = params.bc + params.umi
difference = line.length() - fulllen

}

if ( difference > 0 ) {
println difference
}

process echo {
input:
val difference

when:
difference > 0

script:
"""
echo "Hi: " $difference
"""

}

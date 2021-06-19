params.reads = "fake{1,2}.fastq"
params.R1_barcode = 'true'
params.bc = 5
params.umi = 5
params.fromSTARouts = 'false'

test_ch = Channel.fromPath( params.reads )

reads =  test_ch.toSortedList().get()

if ( params.R1_barcode == 'true' && params.fromSTARouts == 'false' ) {

myFile = new File("${reads[0]}")
String line = myFile.readLines().get(1)
fulllen = params.bc + params.umi
difference = line.length() - fulllen

}

else if ( params.R1_barcode != 'true' && params.fromSTARouts == 'false') {

myFile = new File("${reads[1]}")
String line = myFile.readLines().get(1)
fulllen = params.bc + params.umi
difference = line.length() - fulllen

}

if ( difference > 0  && params.fromSTARouts == 'false') {
println difference
}

process echo {
input:
val difference
file(read) from reads

output:
"R*_*.gz" into read_things_ch

when:
difference > 0 && params.fromSTARouts == 'false'

script:
"""
cp ${reads[0]} R1_${reads[0]}
cp ${reads[1]} R2_${reads[0]}
"""

}

if ( difference < 0  && params.fromSTARouts == 'false') {
read_things_ch = Channel.from(reads)
}

println read_things_ch.toSortedList().get()

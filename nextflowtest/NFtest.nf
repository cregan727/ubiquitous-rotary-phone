params.input = 'true'
params.word = 'Claire'

if ( params.input == 'false' ){
word_ch = Channel.empty()

}
else {
word_ch = Channel.value( params.word )
}



process TRUE {

input:

val word from word_ch

output:
stdout to ch

when:
params.input == 'true'

script:
"""
echo ${word}
"""

}


process FALSE {

input:

output:
stdout to chF

when:
params.input == 'false'

script:
"""
echo 'false'
"""

}




ch
.mix(chF)
.subscribe onNext: { println it }, onComplete: { println 'Done' }


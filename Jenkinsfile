#!groovy

pipeline {

    agent {
        label 'docker-host'
    }

    options {
        timeout(time: 1, unit: 'HOURS')
        disableConcurrentBuilds()
    }

    stages {
        stage('Submodules') {
            agent {
                dockerfile {
                    dir '.devcontainer'
                }
            }
            steps {
                sh 'git submodule update --init --recursive'
            }
        }
        stage('Static Code Checking') {
            agent {
                dockerfile {
                    dir '.devcontainer'
                }
            }
            steps {
                sh 'find src -name \"*.py\" | xargs pylint -f parseable | tee pylint.log'
            }
            post {
                always {
                    recordIssues(
                        tool: pyLint(pattern: 'pylint.log'),
                        unstableTotalHigh: 100
                    )
                }
            }
        }
        stage('Unit tests') {
            agent {
                dockerfile {
                    dir '.devcontainer'
                }
            }
            steps {
                sh 'PYTHONPATH=$WORKSPACE/src:$WORKSPACE/libs/sloscillations pytest --junitxml=pyunit.xml'
            }
            post {
                always {
                    recordIssues (
                        tool: junitParser(pattern: "pyunit.xml")
                    )
                }
            }
        }
    }
    post {
        failure {
            mail to: 'bernd.doser@h-its.org', subject: "FAILURE: ${currentBuild.fullDisplayName}", body: "Failed."
        }
    }
}

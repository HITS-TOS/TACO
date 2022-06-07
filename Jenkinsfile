#!groovy

pipeline {

  agent {
    label 'docker-host'
  }

  options {
    timeout(time: 1, unit: 'HOURS')
  }

  stages {
    stage('Unit tests') {
      agent {
        dockerfile {
          dir '.devcontainer'
        }
      }
      steps {
        sh 'PYTHONPATH=$PYTHONPATH:$pwd/src && pytest'
      }
    }
  }
}

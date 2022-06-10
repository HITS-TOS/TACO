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
        sh 'PYTHONPATH=$WORKSPACE/src pytest'
      }
    }
  }
  post {
    failure {
      mail to: 'bernd.doser@h-its.org', subject: "FAILURE: ${currentBuild.fullDisplayName}", body: "Failed."
    }
  }
}

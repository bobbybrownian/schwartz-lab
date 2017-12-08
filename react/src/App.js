import React, { Component } from 'react';
import { BrowserRouter as Router, Route } from 'react-router-dom';
console.log(Route)
import Research from './containers/Research'
import Methods from './containers/Methods'
import Layout from './containers/Layout'
import Contacts from './containers/Contacts'
import LabHistory from './containers/LabHistory'

class App extends Component {
  constructor(props){
    super(props);
    this.state = {
      selectedLabel : "cheese"
    }
  }



  render () {
    let labelDisplay;
    if(this.state.selectedLabel === "cheese" ){
     let labelDisplay = `<Research />`
    }

    return (

      <div>
        <Router>
          <div>
          <div>
            <Layout />
          </div>
          <div>
            <Route exact path="/" component={Research} />
            <Route exact path="/methods" component={Methods} />
            <Route exact path="/history" component={LabHistory} />
            <Route exact path="/contacts" component={Contacts} />
          </div>
          </div>
        </Router>
      </div>
    )
  }
}

export default App

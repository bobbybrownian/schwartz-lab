import React, { Component } from 'react';
import GoogleMapReact from 'google-map-react';
import CityMapPin from '../components/CityMapPin'

class SimpleMap extends Component {
  constructor(props) {
    super(props)
    this.state = {
		center: {lat: 35, lng: 272},
        zoom:4	
    }
  } 

const GoogleMapConfig = {
  key: "AIzaSyDFkVFlpfTKxJFT8wFVKRpQ9Pdbzh8htsU"
}

  
  render() {
	let mapOptions = {
	  scrollwheel: false,
	  zoomControlOptions: {
		position: 'RIGHT_CENTER',
		style: 'SMALL'
	  },
	  draggable: false,
	  rotateControl: false,
	  scaleControl: false,
	  streetViewControl: false,
	  panControl: false,	  
	};
    return (
      <div>
        <div className="googleDivLayer" style={{height: '450px', width: '1100px', border: 'solid 2px black', left: '-100px', position:'relative'}}>
           <GoogleMapReact
           bootstrapURLKeys={GoogleMapConfig}
           defaultCenter={this.state.center}
           defaultZoom={this.state.zoom}
		   options = {mapOptions}
           >
		   
		  <CityMapPin
             lat={45.4}
             lng={286.6}
             city={'Boston, MA (1983-1991)'}
             members={['Membrane-cytoskeleton interactions and the discovery of integrin signaling']}
             lab = {'Harvard Medical School'}
           />
           <CityMapPin
             lat={36.3}
			 lng={240.65}
			 city={'San Diego, CA (1991-2002)'}
             members={['Integrin signaling and mechanotransduction','Basic mechanisms of integrin signaling and application to cancer','Basic mechanisms of mechanosensing and applications to vascular biology']}
             lab = {'Scripps Research Institute'}
           />
           <CityMapPin
             lat={41.5}
             lng={279.5}
             city={'Charlottesville, VA (2002-2012)'}
             members={['Integrin signaling and mechanotransduction','Basic mechanisms of integrin signaling and application to cancer','Basic mechanisms of mechanosensing and applications to vascular biology and disease']}
             lab = {'University of Virginia'}
           />
           <CityMapPin
             lat={44.5}
             lng={285}
             city={'New Haven, CT (2011-present)'}
             members={['Integrins and mechanotransduction','Basic mechanisms of mechanosensing and applications to vascular biology and disease']}
			 lab = {'Yale University'}
           />
           </GoogleMapReact>
         </div>
       </div>
    )
  }
}

export default SimpleMap

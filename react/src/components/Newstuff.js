import React from 'react'

const Newstuff = (props) => {
  // let methodLink = "/methods/" + props.selectedMethod + "/" + props.pdf_url
  //let newstuffPath = "/newstuffs/" + props.img_url
  return (
    <div className='small-12 columns newstuff-section'>
      <p>{props.description}</p>
    </div>
  )
}

export default Newstuff

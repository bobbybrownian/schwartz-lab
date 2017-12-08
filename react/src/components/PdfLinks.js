import React from 'react';
import { Link } from 'react-router-dom';

const PdfLinks = (props) => {

  // debugger;
  let methodLink = "/methods/" + props.selectedMethod + "/" + props.pdf_url

  return (
    <li>
      <a href={methodLink} download>{props.title}</a>
    </li>
  )
}

export default PdfLinks

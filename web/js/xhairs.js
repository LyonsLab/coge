//Global vars
var can_elemt;
var can_context;
var readout_elemt;
var syn_image;
var location_list = [];

function loadstuff(imageURL) {

    //console.log("DEBUG MESSAGE: onLoad worked!");

    console.log(imageURL);
    //Create image object and load the synteny image
    syn_image = new Image();
    syn_image.src = imageURL;
    document.body.appendChild(syn_image);   //Add image to the HTML DOM so the browser gets details on it

    //When the image is finally loaded...
    syn_image.onload = function () {
        //Remove it from the HTML DOM - we're going to redraw it on the canvas anyway.
        document.body.removeChild(syn_image);

        //Get the Canvas element.
        can_elemt = document.getElementById('myCanvas');
        if (!can_elemt || !can_elemt.getContext) {
            alert("Can't find Canvas element in HTML document - this webpage may not work as intended!");
        }

        //Resize the Canvas element as needed to fit the image
        can_elemt.width = syn_image.width;
        can_elemt.height = syn_image.height;

        // Get the Canvas 2d context.
        can_context = can_elemt.getContext('2d');
        if (!can_context) {
            alert("Can't get 2D context for Canvas element - this webpage may not work as intended!");
        }

        //Get the readout element
/*      readout_elemt = document.getElementById('ReadOut');
        if (!readout_elemt) {
            alert("Can't find ReadOut element in HTML document - this webpage may not work as intended!");
        }*/

        //Draw the synteny image to the Canvas so the user has something to look at initally
        can_context.drawImage(syn_image, 0, 0);
    }
}

function trackpointer(local_event) {
    //var cross_color = '#0000ff';
    var cross_color;

    //Clear the screen
    can_context.fillStyle = '#ffffff';
    can_context.fillRect(0, 0, can_elemt.width, can_elemt.height);
    can_context.drawImage(syn_image, 0, 0);

    //Get offset info from jQuery - the usual method of checking the element offset will not work
    var jq_offset = jQuery('#myCanvas').offset();
    var mp_x = local_event.clientX - jq_offset.left + window.pageXOffset;
    var mp_y = local_event.clientY - jq_offset.top + window.pageYOffset;

    var list_output = list_hitdetect(mp_x, mp_y);
    //If there is no mouseover detail provided, don't color the xhairs. Otherwise, make then red.
    if (list_output[1] == undefined || list_output[1] == '')
    {
        cross_color = '#0000ff';        //Reset to normal blue
    }
    else
    {
        cross_color = '#ff0000';        //Update the crosshair color to red
    }

    //Set color and line width
    can_context.strokeStyle = cross_color;
    can_context.lineWidth   = 0.5;

    //Draw the vertical line
    can_context.beginPath();
    can_context.moveTo(mp_x,0);
    can_context.lineTo(mp_x,can_elemt.height);
    can_context.closePath();
    can_context.stroke();

    //Draw the horizontal line
    can_context.beginPath();
    can_context.moveTo(0, mp_y);
    can_context.lineTo(can_elemt.width, mp_y);
    can_context.closePath();
    can_context.stroke();
}

function getHit(local_event) {
    var jq_offset = jQuery('#myCanvas').offset();

    //ClientWhatever contains the poistion of the cursor in the window
    //offsetWhatever contains the upper left corner of the object in question - in this case, the canvas.
    //Pageoffset accounts for windows scroll.

    //ClientWhatever contains the poistion of the cursor in the window
    //Pageoffset accounts for windows scroll.

    var mp_x = local_event.clientX - jq_offset.left + window.pageXOffset;
    var mp_y = local_event.clientY - jq_offset.top + window.pageYOffset;
    //console.log(local_event, can_elemt.offsetTop);

    //Whenever the mouse is moved, get the mouseOver Javascript and run it
    var list_output = list_hitdetect(mp_x, mp_y);

    //Make sure the list isn't empty
    if (list_output.length) {
        location.href="javascript:"+list_output[1];   //Just want the link, not the mouseover.
    }
}

function trackclick(local_event, new_window) {

    //Get offset info from jQuery - the usual method of checking the element offset will not work
    var jq_offset = jQuery('#myCanvas').offset();
//  console.log(jq_offset); //DEBUG

    //ClientWhatever contains the poistion of the cursor in the window
    //ERIC, Brent added pageOffset here
    var mp_x = local_event.clientX - jq_offset.left + window.pageXOffset;
    var mp_y = local_event.clientY - jq_offset.top + window.pageYOffset;

    //  console.log("User clicked on: "+mp_x+", "+mp_y);
    //alert("Clicked: '"+list_hitdetect(mp_x, mp_y)+"'");

    //Get whatever HREF data is in the massive list and run it.
    var list_output = list_hitdetect(mp_x, mp_y);
//  console.log(hitdetect);     //DEBUG
    if (list_output != [] && list_output[0] != undefined)   //Make sure the list isn't empty
    {
        if (new_window == "no")
        {
            //Open the link in the local window
            location.href=list_output[0];   //Just want the link, not the mouseover.
        }
        else if (new_window == "yes")
        {
            //Open a new window and run the HREF there
            window.open(list_output[0])
        }
    }
}

//Given an x,y position, this function searches the entire list of interesting locations for a point that is within range.
function list_hitdetect(mp_x, mp_y)
{
    var array_length = location_list.length;
    var x_dist;
    var y_dist;

    for (var index = 0; index < array_length; index++)
    {
        if (location_list[index][0] == "circle")
        {
            //Find the distance between the mouse pointer and the current element at the array index
            x_dist = Math.abs(location_list[index][1][0] - mp_x);
            y_dist = Math.abs(location_list[index][1][1] - mp_y);
            var dist = Math.sqrt(Math.pow(x_dist,2)+Math.pow(y_dist,2));

            //DEBUG:
            //readout_elemt.innerHTML = "X:"+mp_x+", Y:"+mp_y+"; Dist:"+dist;

            //If the distance between the mouse pointer and the PoI is less than or equal to the size of the PoI, return the associated info.
            if (dist <= location_list[index][1][2])
            {
                return [location_list[index][2], location_list[index][3]];  //HREF, and OnMouseOver args
            }
        }
        else if (location_list[index][0] == "rect")
        {
            //This is a bit simpler - just a bunch of if statements.
            //NOTE: It does assume that point 1 is closer to the origin than point 2.

            if (mp_x > location_list[index][1][0] && mp_x < location_list[index][1][2])
            {
                if (mp_y > location_list[index][1][1] && mp_y < location_list[index][1][3])
                {
                    //                                      console.log(mp_x, mp_y, location_list[index][1]);
                    //It appears to be inside the square the points make. Go ahead and return the needed info.
                    return [location_list[index][2], location_list[index][3]];  //HREF, and OnMouseOver args
                }
            }
        }
        else
        {
            alert("Unknown areamap type '"+location_list[index][0]+"' at list index "+index);
        }
    }

    //If we get here, there were no hits. Ergo, return an empty list.
    return [];
}
